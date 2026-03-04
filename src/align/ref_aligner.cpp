#include "align.h"
#include "config.hpp"
#include "preprocess.h"
#include "consensus.h"
#include "seed.h"
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include <unordered_map>
#include <cstdio>

namespace align {

    // ==================================================================
    // 构造函数：直接参数初始化
    // 功能：读取参考序列，计算 sketch/minimizer 索引，生成共识序列
    // ==================================================================
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                           int kmer_size, int window_size,
                           int sketch_size, bool noncanonical,
                           int threads, std::string msa_cmd,
                           bool keep_length)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_cmd(std::move(msa_cmd)),
          keep_length(keep_length)
    {
        // 1. 读取参考序列并构建索引
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        while (reader.next(rec))
        {
            // 注意：必须在 move(rec) 之前计算 sketch/minimizer
            auto sketch = mash::sketchFromSequence(rec.seq, kmer_size, sketch_size,
                                                   noncanonical, random_seed);
            auto minimizer = minimizer::extractMinimizer(rec.seq, kmer_size,
                                                         window_size, noncanonical);

            // 保存参考序列及其索引（同索引对应）
            ref_sequences.push_back(std::move(rec));
            ref_sketch.push_back(std::move(sketch));
            ref_minimizers.push_back(std::move(minimizer));
        }

        // 2. 准备共识序列生成的文件路径
        const FilePath consensus_unaligned_file = ref_fasta_path;
        const FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
        const FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

        // 3. 执行 MSA 并生成共识序列
        constexpr std::size_t consensus_batch_size = 4096;
        alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, this->msa_cmd, threads);
        std::string consensus_string = consensus::generateConsensusSequence(
            consensus_aligned_file, consensus_file, consensus_json_file,
            0, threads, consensus_batch_size  // 0 表示使用全部序列
        );

        // 4. 保存共识序列
        consensus_seq.id = "consensus";
        consensus_seq.seq = std::move(consensus_string);

        // 5. 预计算共识序列的 sketch 和 minimizer（性能优化）
        // 用途：避免在每次 alignOneQueryToRef 中重复计算
        consensus_sketch = mash::sketchFromSequence(
            consensus_seq.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        consensus_minimizer = minimizer::extractMinimizer(
            consensus_seq.seq,
            kmer_size,
            window_size,
            noncanonical);
    }

    // ==================================================================
    // 构造函数：基于 Options 初始化
    // 功能：从配置对象中提取参数，委托给主构造函数
    // ==================================================================
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,
            ref_fasta_path,
            opt.kmer_size,
            opt.kmer_window,
            opt.sketch_size,
            true,  // noncanonical 固定为 true
            opt.threads,
            opt.msa_cmd,
            opt.keep_length)
    {
        // 委托构造函数已完成所有初始化
    }

    // ==================================================================
    // 全局序列比对
    // 功能：根据相似度自动选择比对算法
    // ==================================================================
    cigar::Cigar_t RefAligner::globalAlign(const std::string& ref,
                                           const std::string& query,
                                           double similarity,
                                           const SeedHits* ref_minimizer,
                                           const SeedHits* query_minimizer) const
    {
        // 1. 确保 minimizer 可用（如果为空则现场计算）
        const SeedHits* ref_mz_ptr = ref_minimizer;
        const SeedHits* qry_mz_ptr = query_minimizer;

        SeedHits ref_mz_tmp;
        SeedHits qry_mz_tmp;

        if (ref_mz_ptr == nullptr || ref_mz_ptr->empty()) {
            ref_mz_tmp = minimizer::extractMinimizer(ref, kmer_size, window_size, noncanonical);
            ref_mz_ptr = &ref_mz_tmp;
        }
        if (qry_mz_ptr == nullptr || qry_mz_ptr->empty()) {
            qry_mz_tmp = minimizer::extractMinimizer(query, kmer_size, window_size, noncanonical);
            qry_mz_ptr = &qry_mz_tmp;
        }

        // 2. 生成锚点并执行比对
        const anchor::Anchors anchors = minimizer::collect_anchors(*ref_mz_ptr, *qry_mz_ptr);
        cigar::Cigar_t result = globalAlignMM2(ref, query, anchors);

#ifdef _DEBUG
        // Debug 校验：CIGAR 长度必须与序列长度一致
        const std::size_t cigar_ref_len = cigar::getRefLength(result);
        const std::size_t cigar_qry_len = cigar::getQueryLength(result);
        if (cigar_ref_len != ref.size() || cigar_qry_len != query.size()) {
            spdlog::debug("globalAlign: CIGAR 长度不匹配!");
            spdlog::debug("  ref: {}, CIGAR ref: {}", ref.size(), cigar_ref_len);
            spdlog::debug("  query: {}, CIGAR query: {}", query.size(), cigar_qry_len);
            spdlog::debug("  anchors: {}", anchors.size());
        }
#endif
        return result;
    }

    // ==================================================================
    // 写入 SAM 记录
    // ==================================================================
    void RefAligner::writeSamRecord(const seq_io::SeqRecord& q,
                                    const cigar::Cigar_t& cigar,
                                    std::string_view ref_name,
                                    seq_io::SeqWriter& out) const
    {
        const std::string cigar_str = cigar::cigarToString(cigar);
        const auto sam_rec = seq_io::makeSamRecord(q, ref_name, cigar_str, 1, 60, 0);
        out.writeSam(sam_rec);
    }

    // ==================================================================
    // 合并共识序列和 SAM 文件为 FASTA
    // 功能：将共识序列和多个 SAM 文件合并为一个 FASTA 文件
    // 性能：使用大缓冲区批量写入，流式处理避免内存峰值
    // ==================================================================
    std::size_t RefAligner::mergeConsensusAndSamToFasta(
        const std::vector<FilePath>& sam_paths,
        const FilePath& fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map,
        bool keep,
        std::size_t line_width) const
    {
        // 1. 创建 FASTA writer 并写入共识序列
        seq_io::SeqWriter writer(fasta_path, line_width);
        writer.writeFasta(consensus_seq);
        writer.flush();

        std::size_t total_count = 1;  // 已写入共识序列
        std::size_t file_idx = 0;

        // 2. 逐个处理 SAM 文件
        for (const auto& sam_path : sam_paths) {
            // 检查文件是否存在且非空
            if (!std::filesystem::exists(sam_path)) {
                const std::string err_msg = "SAM file does not exist: " + sam_path.string();
                spdlog::error(err_msg);
                throw std::runtime_error(err_msg);
            }

            const auto file_size = std::filesystem::file_size(sam_path);
            if (file_size == 0) {
                spdlog::warn("SAM 文件为空，跳过: {}", sam_path.string());
                ++file_idx;
                continue;
            }

            // 打开并处理当前 SAM 文件
            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;
            seq_io::SeqRecord fasta_rec;
            std::size_t file_count = 0;

            while (reader.next(sam_rec)) {
                // SAM 转 FASTA
                fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);

                // 如果需要保留长度，根据 CIGAR 调整序列
                if (keep && !sam_rec.cigar.empty()) {
#ifdef _DEBUG
                    // Debug：检查 CIGAR 长度一致性
                    const cigar::Cigar_t debug_cigar_ops = cigar::stringToCigar(sam_rec.cigar);
                    const std::size_t debug_ref_len = cigar::getRefLength(debug_cigar_ops);
                    const std::size_t debug_qry_len = cigar::getQueryLength(debug_cigar_ops);
                    if (debug_ref_len != consensus_seq.seq.size() ||
                        debug_qry_len != fasta_rec.seq.size()) {
                        spdlog::debug("mergeConsensusAndSamToFasta: CIGAR 长度不匹配 {}",
                                     fasta_rec.id);
                    }
#endif
                    // 删除 query 相对 ref 的插入，使其投影到参考坐标
                    cigar::Cigar_t cigar_ops = cigar::stringToCigar(sam_rec.cigar);
                    cigar::delQueryToRefByCigar(fasta_rec.seq, cigar_ops);

                    // 如果比对的不是共识序列，需要进一步调整
                    if (sam_rec.rname != "consensus") {
                        auto it = ref_aligned_map.find(sam_rec.rname);
                        if (it == ref_aligned_map.end()) {
                            throw std::runtime_error(
                                "参考序列 '" + sam_rec.rname + "' 未找到");
                        }
                        cigar::padQueryToRefByCigar(fasta_rec.seq, it->second);
                    }

#ifdef _DEBUG
                    // Debug：投影后长度应与共识序列一致
                    if (fasta_rec.seq.size() != consensus_seq.seq.size()) {
                        spdlog::debug("mergeConsensusAndSamToFasta: 投影后长度不匹配 {}",
                                     fasta_rec.id);
                    }
#endif
                }

                // 写入 FASTA
                writer.writeFasta(fasta_rec);
                ++file_count;
                ++total_count;
            }

            ++file_idx;
        }

        // 3. 确保数据写入磁盘
        writer.flush();
        return total_count;
    }

    // ==================================================================
    // 单条 query 比对
    // 功能：计算 query 与最相似参考序列的比对，根据插入情况决定输出
    // ==================================================================
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
    {
        // 1. 计算 query 的 sketch 和 minimizer
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            q.seq, kmer_size, window_size, noncanonical);

        // 2. 选择最相似的参考序列
        seq_io::SeqRecord best_ref;
        double best_jaccard = -1.0;
        std::size_t best_ref_idx = 0;

        if (ref_sequences.size() > 1) {
            // 多参考序列：线性扫描找最相似的
            for (std::size_t r = 0; r < ref_sketch.size(); ++r) {
                const double j = mash::jaccard(qsk, ref_sketch[r]);
                if (j > best_jaccard) {
                    best_jaccard = j;
                    best_ref_idx = r;
                }
            }
            best_ref = ref_sequences[best_ref_idx];
        } else {
            // 单参考序列：直接使用共识序列
            best_ref = consensus_seq;
            best_jaccard = mash::jaccard(qsk, consensus_sketch);
        }

        // 3. 执行全局比对
        cigar::Cigar_t initial_cigar = globalAlign(
            best_ref.seq,
            q.seq,
            best_jaccard,
            &ref_minimizers[best_ref_idx],
            &query_minimizer);

        // 4. 根据参考序列数量和 keep_length 决定策略
        // 单参考序列：直接使用初始比对结果
        if (ref_sequences.size() == 1) {
            if (cigar::hasInsertion(initial_cigar)) {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out_insertion);
            } else {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out);
            }
            return;
        }

        // 多参考序列 + keep_length：保留与最佳参考的比对结果
        if (keep_length) {
            if (cigar::hasInsertion(initial_cigar)) {
                writeSamRecord(q, initial_cigar, best_ref.id, out_insertion);
            } else {
                writeSamRecord(q, initial_cigar, best_ref.id, out);
            }
            return;
        }

        // 多参考序列 + 非 keep_length：有插入时与共识序列二次比对
        cigar::Cigar_t recheck_cigar;
        if (ref_sequences.size() == 1) {
            recheck_cigar = initial_cigar;
        } else {
            const double consensus_similarity = mash::jaccard(qsk, consensus_sketch);
            recheck_cigar = globalAlign(
                consensus_seq.seq,
                q.seq,
                consensus_similarity,
                &consensus_minimizer,
                &query_minimizer);
        }

        // 根据二次比对结果决定输出
        if (cigar::hasInsertion(recheck_cigar)) {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out);
        }
    }

    // ==================================================================
    // 批量比对 query 序列
    // 功能：并行处理 query 序列，每线程独立输出 SAM 文件
    // 性能：使用 OpenMP 并行 + 批处理 + 每线程独立 writer
    // ==================================================================
    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size)
    {
        // 1. 参数检查
        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: 参考序列为空");
        }

        constexpr std::size_t default_batch_size = 2560;
        if (batch_size == 0) {
            batch_size = default_batch_size;
        }

        // 2. 设置线程数
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        const FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        spdlog::info("开始比对：{} 线程，批处理大小 {}", nthreads, batch_size);

        // 3. 为每个线程创建独立的输出文件
        outs_path.clear();
        outs_path.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion_path.clear();
        outs_with_insertion_path.resize(static_cast<std::size_t>(nthreads));

        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs;
        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs_with_insertion;
        outs.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion.resize(static_cast<std::size_t>(nthreads));

        for (int tid = 0; tid < nthreads; ++tid) {
            const FilePath out_path = result_dir /
                (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_SAM_SUFFIX);
            const FilePath out_path_insertion = result_dir /
                (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_INSERTION_SAM_SUFFIX);

            outs_path[static_cast<std::size_t>(tid)] = out_path;
            outs_with_insertion_path[static_cast<std::size_t>(tid)] = out_path_insertion;

            auto tmp = seq_io::SeqWriter::Sam(out_path);
            auto tmp_insertion = seq_io::SeqWriter::Sam(out_path_insertion);

            outs[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp));
            outs[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");

            outs_with_insertion[static_cast<std::size_t>(tid)] =
                std::make_unique<seq_io::SeqWriter>(std::move(tmp_insertion));
            outs_with_insertion[static_cast<std::size_t>(tid)]->writeSamHeader("@HD\tVN:1.6\tSO:unknown");
        }

        // 4. 流式读取 + 批处理并行
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        chunk.reserve(batch_size);

        ProgressBar progress("align");

        while (true) {
            // 释放上一批次的内存（避免峰值内存持续占用）
            chunk.clear();
            chunk.shrink_to_fit();
            chunk.reserve(batch_size);

            // 读取一个批次
            seq_io::SeqRecord rec;
            for (std::size_t i = 0; i < batch_size; ++i) {
                if (!reader.next(rec)) break;
                chunk.push_back(std::move(rec));
            }
            if (chunk.empty()) break;

            // 并行处理当前批次
            #pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk)
            {
                const int tid = omp_get_thread_num();
                auto& out = *outs[static_cast<std::size_t>(tid)];
                auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];

                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i) {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], out, out_insertion);
                }
            }

            // 保存批次大小（用于进度统计）
            const std::size_t chunk_size = chunk.size();

            // 刷新所有 writer
            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }

            // 立即释放内存
            std::vector<seq_io::SeqRecord>().swap(chunk);

            // 更新进度
            progress.tick(chunk_size);
        }

        // 5. 完成并确保所有数据写入磁盘
        progress.done();
        spdlog::info("比对完成");

        for (auto& w : outs) {
            if (w) w->flush();
        }
        for (auto& w : outs_with_insertion) {
            if (w) w->flush();
        }
    }

    // ==================================================================
    // 解析对齐后的参考序列为 CIGAR
    // 功能：读取 MSA 结果，为每条序列生成 CIGAR（M/D 操作）
    // ==================================================================
    void RefAligner::parseAlignedReferencesToCigar(
        const FilePath& aligned_fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
        std::vector<bool>& out_ref_gap_pos) const
    {
        out_ref_aligned_map.clear();
        out_ref_gap_pos.clear();

        seq_io::KseqReader reader(aligned_fasta_path);
        seq_io::SeqRecord rec;
        std::size_t parsed_seq_count = 0;

        while (reader.next(rec)) {
            // 第一条序列：记录 gap 位置
            if (parsed_seq_count == 0) {
                out_ref_gap_pos.reserve(rec.seq.size());
                for (const char base : rec.seq) {
                    out_ref_gap_pos.push_back(base == '-');
                }
            }

            // 为当前序列生成 CIGAR（使用游程编码）
            cigar::Cigar_t cigar;
            cigar.reserve(20);

            char current_op = '\0';
            std::uint32_t current_len = 0;

            for (const char base : rec.seq) {
                const char op = (base == '-') ? 'D' : 'M';
                if (op == current_op) {
                    ++current_len;
                } else {
                    if (current_op != '\0' && current_len > 0) {
                        cigar.push_back(cigar::cigarToInt(current_op, current_len));
                    }
                    current_op = op;
                    current_len = 1;
                }
            }

            // 收尾：写入最后一段
            if (current_op != '\0' && current_len > 0) {
                cigar.push_back(cigar::cigarToInt(current_op, current_len));
            }

            out_ref_aligned_map[rec.id] = std::move(cigar);
            ++parsed_seq_count;
        }

        // 至少要有一条序列
        if (parsed_seq_count == 0) {
            throw std::runtime_error(
                "parseAlignedReferencesToCigar: 输入 FASTA 为空: " +
                aligned_fasta_path.string());
        }

#ifdef _DEBUG
        spdlog::info("parseAlignedReferencesToCigar: 解析 {} 条序列，gap 位置数组长度 {}",
                    parsed_seq_count, out_ref_gap_pos.size());
#endif
    }

    // ==================================================================
    // 辅助函数：将 SAM 记录转换为 FASTA 并根据 CIGAR 调整
    // 功能：单条 SAM 转 FASTA 的转换逻辑，支持 CIGAR 对齐
    // ==================================================================
    void RefAligner::convertSamToFastaRecord(
        const seq_io::SamRecord& sam_rec,
        seq_io::SeqRecord& fasta_rec,
        const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
        std::size_t estimated_final_length) const
    {
        // 转换 SAM 为 FASTA
        fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);

        // 预分配内存以提高性能
        if (fasta_rec.seq.capacity() < estimated_final_length) {
            fasta_rec.seq.reserve(estimated_final_length);
        }

        // 根据 SAM 的 CIGAR 字段调整序列
        if (!sam_rec.cigar.empty() && sam_rec.cigar != "*") {
            cigar::Cigar_t tmp_cigar = cigar::stringToCigar(sam_rec.cigar);
            cigar::padQueryToRefByCigar(fasta_rec.seq, tmp_cigar);
        }

        // 根据参考序列的对齐信息进一步调整
        auto it = ref_aligned_map.find(sam_rec.rname);
        if (it == ref_aligned_map.end()) {
            throw std::runtime_error(
                "Reference sequence '" + sam_rec.rname + "' not found in alignment map");
        }
        cigar::padQueryToRefByCigar(fasta_rec.seq, it->second);
    }

    // ==================================================================
    // 辅助函数：处理插入序列文件
    // 功能：读取含插入的 SAM，合并为 FASTA，执行可选 MSA
    // ==================================================================
    FilePath RefAligner::processInsertionSequences(
        const FilePath& result_dir,
        const FilePath& aligned_insertion_fasta,
        std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map) const
    {
        // 收集插入 SAM 文件路径
        std::vector<FilePath> insertion_sam_paths(outs_with_insertion_path.begin(),
                                                   outs_with_insertion_path.end());
        spdlog::info("Collecting insertion SAM files: {} files", insertion_sam_paths.size());

        const FilePath insertion_fasta_path = result_dir / ALL_INSERTION_FASTA;
        const std::size_t total_sequences = mergeConsensusAndSamToFasta(
            insertion_sam_paths,
            insertion_fasta_path,
            ref_aligned_map,
            keep_length,
            80);

        spdlog::info("Merged insertion sequences: {} records (including consensus)", total_sequences);

        // 执行 MSA 或复制文件
        if (!keep_length) {
            spdlog::info("Running MSA on insertion sequences");
            alignConsensusSequence(insertion_fasta_path, aligned_insertion_fasta,
                                  msa_cmd, threads);
        } else {
            file_io::copyFile(insertion_fasta_path, aligned_insertion_fasta);
            spdlog::info("Skipped MSA (keep_length=true), file copied");
        }

        return insertion_fasta_path;
    }

    // ==================================================================
    // 辅助函数：写入共识和参考序列
    // 功能：从对齐文件读取并写入共识及参考序列
    // ==================================================================
    std::size_t RefAligner::writeConsensusAndReferences(
        seq_io::SeqWriter& final_writer,
        const FilePath& consensus_aligned_file,
        ProgressBar& progress) const
    {
        spdlog::info("Writing consensus and reference sequences");
        seq_io::KseqReader cons_reader(consensus_aligned_file);
        seq_io::SeqRecord cons_rec;
        std::size_t seq_count = 0;

        while (cons_reader.next(cons_rec)) {
            seq_io::cleanSequence(cons_rec);
            final_writer.writeFasta(cons_rec);
            ++seq_count;
            progress.tick();
        }
        final_writer.flush();

        spdlog::info("Consensus and reference sequences written: {} records", seq_count);
        return seq_count;
    }

    // ==================================================================
    // 辅助函数：写入插入序列
    // 功能：从对齐的插入文件读取并写入序列（跳过第一条共识）
    // ==================================================================
    std::size_t RefAligner::writeInsertionSequences(
        seq_io::SeqWriter& final_writer,
        const FilePath& aligned_insertion_fasta,
        std::size_t& expected_length,
        bool& length_initialized,
        ProgressBar& progress) const
    {
        spdlog::info("Writing insertion sequences");
        seq_io::KseqReader insertion_reader(aligned_insertion_fasta);
        seq_io::SeqRecord insertion_rec;
        std::size_t seq_count = 0;
        bool skip_first = true;

        while (insertion_reader.next(insertion_rec)) {
            // 跳过共识序列（第一条）
            if (skip_first) {
                skip_first = false;
                continue;
            }

            // 检查长度一致性
            if (!length_initialized) {
                expected_length = insertion_rec.seq.size();
                length_initialized = true;
            } else if (insertion_rec.seq.size() != expected_length) {
                throw std::runtime_error(
                    "Sequence length mismatch: " + insertion_rec.id +
                    " length " + std::to_string(insertion_rec.seq.size()) +
                    ", expected " + std::to_string(expected_length));
            }

            final_writer.writeFasta(insertion_rec);
            ++seq_count;
            progress.tick();
        }
        final_writer.flush();

        spdlog::info("Insertion sequences written: {} records total", seq_count);
        return seq_count;
    }

    // ==================================================================
    // 辅助函数：处理单个 SAM 文件的批次
    // 功能：读取 SAM 批次，并行转换为 FASTA，串行写入输出
    // ==================================================================
    void RefAligner::processSamFileBatch(
        seq_io::SamReader& sam_reader,
        const std::size_t batch_size,
        seq_io::SeqWriter& final_writer,
        const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
        std::size_t estimated_final_length,
        std::size_t& expected_length,
        bool& length_initialized,
        std::size_t& seq_count,
        ProgressBar& progress) const
    {
        std::vector<seq_io::SamRecord> sam_batch;
        std::vector<seq_io::SeqRecord> fasta_batch;
        sam_batch.reserve(batch_size);
        fasta_batch.resize(batch_size);

        seq_io::SamRecord sam_rec;
        while (true) {
            // 读取一个批次
            sam_batch.clear();
            while (sam_batch.size() < batch_size && sam_reader.next(sam_rec)) {
                sam_batch.push_back(sam_rec);
            }

            if (sam_batch.empty()) break;

            const std::size_t current_batch_size = sam_batch.size();

            // 并行转换 SAM 为 FASTA
            #pragma omp parallel for default(none) \
                shared(sam_batch, fasta_batch, current_batch_size, ref_aligned_map, \
                       estimated_final_length) \
                schedule(dynamic, 4) num_threads(threads)
            for (std::size_t i = 0; i < current_batch_size; ++i) {
                convertSamToFastaRecord(sam_batch[i], fasta_batch[i],
                                       ref_aligned_map, estimated_final_length);
            }

            // 串行写入（确保顺序和长度检查）
            for (std::size_t i = 0; i < current_batch_size; ++i) {
                const seq_io::SeqRecord& fasta_rec = fasta_batch[i];

                // 检查长度一致性
                if (!length_initialized) {
                    expected_length = fasta_rec.seq.size();
                    length_initialized = true;
                } else if (fasta_rec.seq.size() != expected_length) {
                    throw std::runtime_error(
                        "Sequence length mismatch: " + fasta_rec.id +
                        " length " + std::to_string(fasta_rec.seq.size()) +
                        ", expected " + std::to_string(expected_length));
                }

                final_writer.writeFasta(fasta_rec);
                ++seq_count;
                progress.tick();
            }
        }
    }

    // ==================================================================
    // 合并比对结果
    // 功能：合并所有线程输出的 SAM 文件为最终 FASTA
    // ==================================================================
    void RefAligner::mergeAlignedResults(const FilePath output, std::size_t batch_size)
    {
        ProgressBar progress("merge");

        const FilePath result_dir = work_dir / RESULTS_DIR;
        const FilePath aligned_insertion_fasta = result_dir / ALIGNED_INSERTION_FASTA;
        const FilePath consensus_aligned_file =
            FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;

        // 1. 解析参考序列的对齐信息
        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map;
        std::unordered_map<std::string, cigar::Cigar_t> insertion_aligned_map;
        std::vector<bool> ref_gap_pos;
        std::vector<bool> insertion_ref_gap_pos;

        parseAlignedReferencesToCigar(consensus_aligned_file, ref_aligned_map, ref_gap_pos);

        // 2. 创建最终输出 writer
        spdlog::info("Writing final MSA FASTA: {}", output.string());
        seq_io::SeqWriter final_writer(output, U_MAX);

        // 3. 处理插入序列
        processInsertionSequences(result_dir, aligned_insertion_fasta, ref_aligned_map);

        // 4. 解析插入序列的对齐信息
        parseAlignedReferencesToCigar(aligned_insertion_fasta, insertion_aligned_map,
                                      insertion_ref_gap_pos);
        ref_aligned_map[consensus_seq.id] = insertion_aligned_map[consensus_seq.id];

        spdlog::info("Parsed MSA alignment maps: ref_map={}, insertion_map={}",
                    ref_aligned_map.size(), insertion_aligned_map.size());

        std::size_t expected_length = 0;
        std::size_t seq_count = 0;
        bool length_initialized = false;

        // 5. 写入共识和参考序列
        seq_count = writeConsensusAndReferences(final_writer, consensus_aligned_file, progress);

        // 6. 写入插入序列
        seq_count += writeInsertionSequences(final_writer, aligned_insertion_fasta,
                                            expected_length, length_initialized, progress);

        // 7. 批处理 SAM 文件
        constexpr std::size_t default_batch_size = 2560;
        const std::size_t effective_batch_size = (batch_size > 0) ? batch_size : default_batch_size;
        const std::size_t estimated_final_length = expected_length > 0 ? expected_length : 30000;

        spdlog::info("Processing SAM files: {} files, batch size {}",
                    outs_path.size(), effective_batch_size);

        for (const auto& sam_path : outs_path) {
            seq_io::SamReader sam_reader(sam_path);
            processSamFileBatch(sam_reader, effective_batch_size, final_writer,
                               ref_aligned_map, estimated_final_length,
                               expected_length, length_initialized, seq_count, progress);
        }

        final_writer.flush();

        // 8. 完成
        progress.done();
        spdlog::info("Merge completed: {} sequences total, length {}",
                    seq_count, expected_length);
    }

    // ==================================================================
    // 删除参考序列的 gap 列
    // 功能：原地过滤，使用双指针优化性能
    // ==================================================================
    void RefAligner::removeRefGapColumns(
        std::string& seq,
        const std::vector<bool>& ref_gap_pos)
    {
        // Fast-path：gap 标记为空时无需处理
        if (ref_gap_pos.empty()) {
            return;
        }

#ifdef _DEBUG
        // Debug：检查长度一致性
        if (seq.size() != ref_gap_pos.size()) {
            throw std::runtime_error(
                "removeRefGapColumns: 序列长度不匹配 seq_len=" +
                std::to_string(seq.size()) +
                ", ref_gap_pos_len=" + std::to_string(ref_gap_pos.size()));
        }
#endif

        // 原地过滤：使用双指针
        std::size_t write_pos = 0;
        const std::size_t n = seq.size();

        for (std::size_t read_pos = 0; read_pos < n; ++read_pos) {
            if (!ref_gap_pos[read_pos]) {
                seq[write_pos++] = seq[read_pos];
            }
        }

        // 截断多余部分
        seq.resize(write_pos);
    }

} // namespace align

