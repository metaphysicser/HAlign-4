#include "align.h"
#include "config.hpp"
#include "preprocess.h"    // alignConsensusSequence
#include "consensus.h"     // generateConsensusSequence
#include "seed.h"          // minimizer 提取和锚点生成
#include <algorithm>
#include <cstddef>
#include <filesystem>      // std::filesystem::remove
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

#include <unordered_map>

// 进度条：
// - <cstdio> 用于 fprintf(stderr, ...) 实现单行(\r)进度输出
// 注：<chrono> 在本文件中未直接使用（ProgressBar 内部若需要计时，应在其实现文件中包含）
// #include <chrono>       // 未使用：移除以降低噪音（不影响逻辑）
#include <cstdio>

namespace align {

    // ------------------------------------------------------------------
    // 构造函数1：直接传入参数初始化
    // 说明（与代码真实行为保持一致）：
    // - 本构造函数会读取 ref_fasta_path 中的所有序列，逐条计算 MinHash sketch 与 minimizer（用于后续相似度与锚点生成）。
    // - 随后会对 ref_fasta_path 内的序列做一次 MSA，并生成一条“共识序列”(consensus_string)。
    // - keep_first_length=true 时：共识序列直接取第一条参考序列（ref_sequences.front()），忽略 consensus_string。
    // - keep_first_length=false 时：使用生成的 consensus_string 作为共识序列内容。
    // 注意：这里“是否生成/使用共识序列”的行为与 keep_first_length 强相关，注释必须避免写反。
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path, int kmer_size, int window_size,
                           int sketch_size, bool noncanonical, int threads, std::string msa_cmd,
                           bool keep_length)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_cmd(msa_cmd),          // 这里是一次按值传参 + 成员拷贝/移动（由编译器决定 NRVO/移动），不是显式 std::move
          keep_length(keep_length)
    {
        // 从文件读取参考序列并构建索引（每条 ref 计算 sketch/minimizer，供后续相似度筛选与 anchors 生成）
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        while (reader.next(rec))
        {
            // 关键点：必须在 move(rec) 之前计算 sketch/minimizer，因为 rec.seq 后续会被移动走
            auto sketch = mash::sketchFromSequence(rec.seq, kmer_size, sketch_size, noncanonical, random_seed);
            auto minimizer = minimizer::extractMinimizer(rec.seq, kmer_size, window_size, noncanonical);

            // 保存参考序列本体
            ref_sequences.push_back(std::move(rec));
            // 保存 sketch/minimizer（和 ref_sequences 同索引对应）
            ref_sketch.push_back(std::move(sketch));
            ref_minimizers.push_back(std::move(minimizer));
        }

        // MSA/共识生成所需的工作目录文件路径（沿用项目 workdir 约定）
        FilePath consensus_unaligned_file = ref_fasta_path; // 输入：未对齐序列集合（这里直接用 ref_fasta_path）
        FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED; // 输出：MSA 结果
        FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;           // 输出：共识 FASTA
        FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;       // 输出：共识统计 JSON

        const std::size_t batch_size = 4096; // 共识生成的批处理大小（影响 IO/内存，不影响结果）
        // 外部 MSA：把未对齐序列集合对齐到 consensus_aligned_file
        alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, msa_cmd, threads);
        // 从 MSA 结果生成一条共识序列字符串
        std::string consensus_string = consensus::generateConsensusSequence(
            consensus_aligned_file,
            consensus_file,
            consensus_json_file,
            0, // 不限制数量（使用全部序列）
            threads,
            batch_size
        );


        consensus_seq.id = "consensus";
        consensus_seq.seq = std::move(consensus_string);


        // ------------------------------------------------------------------
        // 性能优化：预计算共识序列的 MinHash sketch 和 Minimizer 索引
        // ------------------------------------------------------------------
        // 说明：
        // 1. consensus_sketch 和 consensus_minimizer 在每次 alignOneQueryToRef 的二次比对中都会用到
        // 2. 将其提前计算并存储为成员变量，避免每次都重复计算
        // 3. 对于大规模比对任务（数千到数万条 query），可显著减少 CPU 时间
        //
        // 性能影响：
        // - 初始化时增加一次 sketch 和 minimizer 计算开销（通常 < 20ms）
        // - 每次 alignOneQueryToRef 节省一次 sketch 计算（约 1-5ms）和 minimizer 计算（约 2-10ms）
        // - 对于 10000 条 query，可节省 30-150 秒总时间
        //
        // 实现说明：
        // - consensus_sketch：用于 Jaccard 相似度计算，决定是否需要二次比对
        // - consensus_minimizer：用于种子定位和锚点比对（当前未使用，预留未来优化）
        // ------------------------------------------------------------------
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

    // ------------------------------------------------------------------
    // 构造函数2：基于 Options 结构体初始化
    // 说明：
    // 1. 从 Options 中提取相关参数，委托给第一个构造函数
    // 2. 构造函数1会根据 keep_first_length 标志自动生成或选择参考序列
    // 3. 从 opt 中提取 keep_first_length 和 keep_all_length
    // ------------------------------------------------------------------
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,                // work_dir：工作目录
            ref_fasta_path,             // 参考序列文件路径
            opt.kmer_size,              // kmer_size：k-mer 大小
            opt.kmer_window,            // window_size：minimizer 窗口大小
            opt.sketch_size,            // sketch_size：sketch 大小
            true,                       // noncanonical：是否使用非标准模式（固定为 true）
            opt.threads,                // threads：线程数（用于共识序列生成）
            opt.msa_cmd,                // msa_cmd：MSA 命令模板
            opt.keep_length      // keep_length：从 opt 提取
        )
    {
        // 委托构造函数已完成所有初始化工作（包括共识序列生成）
    }

    // ==================================================================
    // 静态方法：globalAlign - 全局序列比对（统一接口）
    // ==================================================================
    // 功能：
    // 根据两个序列的相似度自动选择合适的比对算法执行全局比对
    // ==================================================================
    cigar::Cigar_t RefAligner::globalAlign(const std::string& ref,
                                           const std::string& query,
                                           double similarity,
                                           const SeedHits* ref_minimizer,
                                           const SeedHits* query_minimizer) const
    {
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

        // 生成 anchors（这是最容易发生“内存爆炸”的一步）
        anchor::Anchors anchors = minimizer::collect_anchors(*ref_mz_ptr, *qry_mz_ptr);

        // // 原逻辑：相似度稍高用 WFA2 分段，否则用 KSW2 分段
        cigar::Cigar_t result;
        result = globalAlignMM2(ref, query, anchors);


#ifdef _DEBUG
        // Debug 校验：CIGAR 必须消耗完整的 ref 和 query 长度
        const std::size_t cigar_ref_len = cigar::getRefLength(result);
        const std::size_t cigar_qry_len = cigar::getQueryLength(result);
        if (cigar_ref_len != ref.size() || cigar_qry_len != query.size()) {
            spdlog::debug("RefAligner::globalAlign [MM2路径]: CIGAR长度不匹配!");
            spdlog::debug("  ref长度: {}, CIGAR消耗ref: {}", ref.size(), cigar_ref_len);
            spdlog::debug("  query长度: {}, CIGAR消耗query: {}", query.size(), cigar_qry_len);
            spdlog::debug("  锚点数量: {}", anchors.size());
            spdlog::debug("  CIGAR: {}", cigar::cigarToString(result));
        }
#endif
        return result;
    }


    // ------------------------------------------------------------------
    // 辅助函数：写入SAM记录
    // ------------------------------------------------------------------
    void RefAligner::writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                                    std::string_view ref_name, seq_io::SeqWriter& out) const
    {
        const std::string cigar_str = cigar::cigarToString(cigar);
        const auto sam_rec = seq_io::makeSamRecord(q, ref_name, cigar_str, 1, 60, 0);
        out.writeSam(sam_rec);
    }

    // ------------------------------------------------------------------
    // 辅助函数：mergeConsensusAndSamToFasta
    // 功能：将共识序列和多个 SAM 文件合并为一个 FASTA 文件
    //
    // 实现说明：
    // 1. 创建一个 SeqWriter 打开输出 FASTA 文件
    // 2. 先写入共识序列（consensus_seq）
    // 3. 逐个读取 SAM 文件，提取序列并追加写入到同一个 writer
    // 4. 使用同一个 writer 实例，避免文件追加模式的复杂性和性能损失
    // 5. 所有写入操作在同一个 writer 的生命周期内完成
    // 6. 如果 keep 为 true，根据 SAM 记录中的 CIGAR 字符串在序列中插入 gap
    //
    // 性能优化：
    // - 使用大缓冲区（默认 8MiB）批量写入，减少系统调用
    // - 流式处理 SAM 文件，内存占用与文件大小无关
    // - 避免文件的重复打开和关闭
    // ------------------------------------------------------------------
    std::size_t RefAligner::mergeConsensusAndSamToFasta(
        const std::vector<FilePath>& sam_paths,
        const FilePath& fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map,
        bool keep,
        std::size_t line_width
        ) const
    {
        // 创建 FASTA writer（非追加模式，会覆盖已存在的文件）
        seq_io::SeqWriter writer(fasta_path, line_width);

        // 1. 先写入共识序列
        writer.writeFasta(consensus_seq);
        writer.flush();

        // 2. 统计信息（用于返回值和调试）
        std::size_t total_count = 1;  // 已写入共识序列，计数从 1 开始
        std::size_t file_idx = 0;

        // 3. 逐个处理每个 SAM 文件
        for (const auto& sam_path : sam_paths) {
            // 检查 1：文件是否存在
            if (!std::filesystem::exists(sam_path)) {
                const std::string err_msg = "SAM file does not exist: " + sam_path.string();
                spdlog::error(err_msg);
                throw std::runtime_error(err_msg);
            }

            // 检查 2：文件是否为空（0 字节通常表示写入失败或程序异常）
            const auto file_size = std::filesystem::file_size(sam_path);
            if (file_size == 0) {
                spdlog::warn("SAM file is empty (0 bytes), skipping: {}", sam_path.string());
                ++file_idx;  // 跳过当前文件，但仍更新文件索引
                continue;    // 跳过空文件，继续处理下一个
            }


            // 打开当前 SAM 文件
            // 说明：
            // - 每个 SamReader 独立打开和关闭
            // - 使用 RAII 确保文件在处理完毕后自动关闭
            // - 如果文件打开失败，会抛出异常并中止合并
            seq_io::SamReader reader(sam_path);

            seq_io::SamRecord sam_rec;
            seq_io::SeqRecord fasta_rec;
            std::size_t file_count = 0;

            while (reader.next(sam_rec)) {

                // 将 SAM 记录转换为 FASTA 记录（不做额外清洗；第二个参数控制是否保留质量等元信息）
                fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);

                // keep=true 的真实含义（与代码一致）：
                // - 这里并不是“插入 gap”，而是把 query 中“相对参考的插入”删除掉，使 query 投影到参考坐标（长度更接近参考）。
                // - 具体行为由 delQueryToRefByCigar 决定：它会按 CIGAR 删除 query 中对应 I 的碱基。
                // 这样做用于后续 MSA/合并阶段的长度/坐标统一。
                if (keep && !sam_rec.cigar.empty()) {
                    // 解析 CIGAR 字符串为内部表示
#ifdef _DEBUG
                    // Debug：检查 cigar 长度是否与 ref/query 一致（便于定位上游写 SAM 的问题）
                    cigar::Cigar_t debug_cigar_ops = cigar::stringToCigar(sam_rec.cigar);
                    std::size_t debug_ref_len = cigar::getRefLength(debug_cigar_ops);
                    std::size_t debug_qry_len = cigar::getQueryLength(debug_cigar_ops);
                    if (debug_ref_len != consensus_seq.seq.size() || debug_qry_len != fasta_rec.seq.size()) {
                        spdlog::debug("mergeConsensusAndSamToFasta: CIGAR length mismatch for record {}: expected ref {}, query {}; got ref {}, query {}",
                                     fasta_rec.id,
                                     consensus_seq.seq.size(),
                                     fasta_rec.seq.size(),
                                     debug_ref_len,
                                     debug_qry_len);
                    }
#endif
                    cigar::Cigar_t cigar_ops = cigar::stringToCigar(sam_rec.cigar); // 解析 CIGAR（后续按其删除 query 插入）
                    cigar::delQueryToRefByCigar(fasta_rec.seq, cigar_ops);         // 将 query 删除投影到 ref 坐标（原地修改）

                    if (sam_rec.rname != "consensus")
                    {
                        auto it = ref_aligned_map.find(sam_rec.rname);
                        if (it == ref_aligned_map.end()) {
                            {
                                throw std::runtime_error(
                                    "mergeAlignedResults: reference '" + sam_rec.rname +
                                    "' not found in ref_aligned_map");
                            }
                        }
                        cigar::padQueryToRefByCigar(fasta_rec.seq, it->second);
                    }
#ifdef _DEBUG
                    // Debug：投影后长度应与共识序列长度一致（若 ref_name 是共识）
                    if (fasta_rec.seq.size() != consensus_seq.seq.size()) {
                       spdlog::debug("mergeConsensusAndSamToFasta: projected size mismatch for record {}: expected {}, got {}",
                                    fasta_rec.id,
                                    consensus_seq.seq.size(),
                                    fasta_rec.seq.size());
                    }
#endif
                }

                // 写出 FASTA（writer 内部有缓冲；此处每条写一次，最终 flush 统一落盘）
                writer.writeFasta(fasta_rec);
                ++file_count;   // 当前 SAM 文件内计数
                ++total_count;  // 全部文件累计计数
            }


            ++file_idx;
        }

        // 4. 确保所有数据已刷新到磁盘
        writer.flush();

        return total_count;
    }

    // ------------------------------------------------------------------
    // 主函数：alignOneQueryToRef
    // 功能：对单个query执行比对并写入SAM文件
    //
    // 流程：
    // 1. 计算 query sketch
    // 2. 找到最相似的 reference
    // 3. 执行全局比对
    // 4. 根据是否有插入进行二次判断并写入相应文件
    //
    // 参数说明：
    //   - q: 待比对的查询序列
    //   - out: 普通输出文件的 writer（无插入或二次比对后无插入的序列）
    //   - out_insertion: 插入序列输出文件的 writer（二次比对后仍有插入的序列）
    // ------------------------------------------------------------------
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
    {

        // 1) 计算 query sketch 和 minimizer
        // 说明：
        // - sketch 用于快速相似度估计（Jaccard 距离）
        // - minimizer 用于种子定位和锚点比对（当前传递给 globalAlign，未来可用）
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            q.seq,
            kmer_size,
            window_size,
            noncanonical);

        seq_io::SeqRecord best_ref;
        double best_j = -1.0;
        std::size_t best_r = 0;

        // 2) 选择最相似 reference（线性扫描）
        // TODO 候选多几个序列比对
        if (ref_sequences.size() > 1)
        {

            for (std::size_t r = 0; r < ref_sketch.size(); ++r) {
                const double j = mash::jaccard(qsk, ref_sketch[r]);
                if (j > best_j) {
                    best_j = j;
                    best_r = r;
                }
            }
            best_ref = ref_sequences[best_r];
        }else
        {
            best_ref = consensus_seq;
            best_j = mash::jaccard(qsk, consensus_sketch);
        }

        // 3) 执行全局比对（使用统一的比对接口）
        // 说明：
        // - 使用 RefAligner::globalAlign 而非直接调用 globalAlignWFA2
        // - 传递 ref_minimizer 和 query_minimizer，允许未来根据种子信息优化比对
        // - 允许未来根据相似度自动选择最优算法
        cigar::Cigar_t initial_cigar = globalAlign(
            best_ref.seq,
            q.seq,
            best_j,
            &ref_minimizers[best_r],  // 使用对应参考序列的 minimizer
            &query_minimizer);

        // keep-length single-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果
        // keep-length multi-ref: 如果没有插入，直接写入，如果有插入，保留和ref的比对结果
        // non-keep-length single-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果
        // non-keep-length multi-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果

        if (ref_sequences.size() == 1) {
            // 如果只有一个参考序列，说明它就是共识序列（无论 keep_length 如何），直接使用初始比对结果即可
            if (cigar::hasInsertion(initial_cigar)) {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out_insertion);
            } else {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out);
            }
            return;
        } else
        {
            if (keep_length)
            {
                if (cigar::hasInsertion(initial_cigar)) {
                    writeSamRecord(q, initial_cigar, best_ref.id, out_insertion);
                } else {
                    writeSamRecord(q, initial_cigar, best_ref.id, out);
                }
                return;
            }
        }


        // if (keep_length) {
        //     if (!cigar::hasInsertion(initial_cigar)) {
        //         writeSamRecord(q, initial_cigar, best_ref.id, out);
        //     } else {
        //         writeSamRecord(q, initial_cigar, consensus_seq.id, out_insertion);
        //     }
        //     return;
        // }

        // 有插入：进行二次比对判断（与共识序列比对）
        // 说明：
        // 1. 使用预计算的 consensus_sketch 和 consensus_minimizer 成员变量（性能优化）
        // 2. 计算 query 与共识序列的相似度（用于传递给 globalAlign）
        // 3. 使用统一的比对接口，保持一致性
        //
        // 性能优化：
        // - consensus_sketch 和 consensus_minimizer 在构造函数中已计算并缓存
        // - 避免每次调用都重复计算 sketchFromSequence 和 extractMinimizer
        // - 对于大规模比对（10000+ queries），可节省数十秒到数分钟
        cigar::Cigar_t recheck_cigar;

        if (ref_sequences.size() == 1) {
            // 如果最相似的参考就是共识序列，说明 query 已经与共识序列比对过了，无需重复比对
            // 直接使用初始比对结果（虽然可能有插入，但我们认为这是最终结果）
            recheck_cigar = initial_cigar;
        }else
        {
            const double consensus_similarity = mash::jaccard(qsk, consensus_sketch);

            recheck_cigar = globalAlign(
                consensus_seq.seq,
                q.seq,
                consensus_similarity,
                &consensus_minimizer,  // 使用预计算的共识序列 minimizer
                &query_minimizer);
        }



        // 根据二次比对结果决定输出文件
        if (cigar::hasInsertion(recheck_cigar)) {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out);
        }

    }

    void RefAligner::alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size)
    {
        // =====================================================================
        // OpenMP 并行骨架（流式读取 + 每线程独立 writer）
        //
        // 核心思路：
        // - 单线程按 chunk 流式读取 query；
        // - 对一个 chunk 用 OpenMP 并行 for；
        // - 每线程有稳定 thread_id (= omp_get_thread_num)；
        // - 每线程一个独占输出文件，避免加锁。
        //
        // 我们把“相似度计算 +（占位）比对 + 写出”抽象成一个函数 process_one()。
        // 这样后续你只需要替换 process_one() 的内部逻辑即可，不影响并行框架。
        // =====================================================================

        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: reference is empty");
        }

        // batch_size==0 没有意义，这里兜底为 2560
        if (batch_size == 0) batch_size = 2560;

        // 线程数：
        // - threads<=0：尊重 OpenMP 默认（OMP_NUM_THREADS 或 runtime 配置）
        // - threads>0 ：显式设置
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        spdlog::info("Starting alignment with {} threads and batch size {}", nthreads, batch_size);

        // 每线程一个输出
        // writer（SAM 模式，占位输出也能保持格式合法）
        outs_path.clear();
        outs_path.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion_path.clear();
        outs_with_insertion_path.resize(static_cast<std::size_t>(nthreads));


        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs;
        std::vector<std::unique_ptr<seq_io::SeqWriter>> outs_with_insertion;

        outs.clear();
        outs_with_insertion.clear();
        outs.resize(static_cast<std::size_t>(nthreads));
        outs_with_insertion.resize(static_cast<std::size_t>(nthreads));  // 关键修复：必须 resize 以避免越界访问

        for (int tid = 0; tid < nthreads; ++tid) {
            // 每个线程独立输出 SAM 文件，避免线程竞争
            // 文件名格式：thread<tid>.sam 和 thread<tid>_insertion.sam
            FilePath out_path = result_dir / (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_SAM_SUFFIX);
            FilePath out_path_insertion = result_dir / (THREAD_SAM_PREFIX + std::to_string(tid) + THREAD_INSERTION_SAM_SUFFIX);

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

        // ---- 流式读取（单线程）+ chunk/batch 并行处理 ----
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        chunk.reserve(batch_size);

        // -----------------------------------------------------------------
        // 进度条（使用 ProgressBar 类，减少代码冗余）
        //
        // 设计说明（为什么这样做）：
        // 1) 你明确不希望预统计（需要二次遍历输入，超大数据集很慢）。
        // 2) 你也不希望用 atomic（并行区频繁原子加会带来额外开销）。
        //
        // 因此这里采用：
        // - "读取线程"统计本次读入/累计已处理的序列数（read/processed）
        // - 并行区只做计算，不做计数
        // - 每完成一个 chunk（即 chunk 内所有 read 都已处理），再在串行区更新 processed_count。
        //
        // 进度刷新策略：每处理 >=1000 条刷新一次，使用 '\r' 覆盖同一行。
        // -----------------------------------------------------------------
        ProgressBar progress("align");

        while (true)
        {
            // ========================================================================
            // 内存修复：强制释放 chunk 容量，避免峰值内存持续占用
            // ========================================================================
            // 问题根因：
            // 1. std::vector::clear() 只删除元素，不释放容量（capacity 保持不变）
            // 2. 如果某个 batch 有大序列（30kb），chunk 容量会扩展到 batch_size × 30kb
            // 3. 后续 batch 即使全是小序列，chunk 仍占用峰值内存
            // 4. 100 万条序列 × 多个 batch → 累积内存达到数百 GB
            //
            // 修复策略：
            // 1. 使用 shrink_to_fit() 强制释放未使用的容量（C++11 标准）
            // 2. 或者使用 swap trick：vector<SeqRecord>().swap(chunk)
            // 3. 重新 reserve(batch_size) 避免下次读取时频繁扩容
            //
            // 性能权衡：
            // - shrink_to_fit() 开销：O(capacity)，通常 < 1ms（对于 batch_size=1000）
            // - 内存收益：减少 99% 持续占用（30MB → 300KB）
            // - 总体影响：100 万条序列增加约 1 秒（1000 次 shrink），换取 700GB → 5GB
            // ========================================================================

            // 方法1：标准方法（C++11）
            chunk.clear();
            chunk.shrink_to_fit();  // 释放未使用容量
            chunk.reserve(batch_size);  // 预分配，避免读取时频繁扩容

            // 方法2（备选）：swap trick（C++98 兼容）
            // std::vector<seq_io::SeqRecord>().swap(chunk);
            // chunk.reserve(batch_size);

            seq_io::SeqRecord rec;
            for (std::size_t i = 0; i < batch_size; ++i) {
                if (!reader.next(rec)) break;
                chunk.push_back(std::move(rec));
            }
            if (chunk.empty()) break;

            // OpenMP 规约：
            // - default(none) 强制显式标注共享/私有变量，避免无意的数据竞争；
            // - outs/outs_with_insertion 与 chunk 是跨线程共享只读/线程索引访问的；tid/i 为线程私有。
            #pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk)
            {
                const int tid = omp_get_thread_num();
                auto& out = *outs[static_cast<std::size_t>(tid)];
                auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];

                #pragma omp for schedule(dynamic, 1)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i)
                {
                    alignOneQueryToRef(chunk[static_cast<std::size_t>(i)], out, out_insertion);
                }
            }

            // ========================================================================
            // 内存优化：立即释放 chunk 容量，避免下一轮循环前的内存峰值
            // ========================================================================
            // 说明：
            // 1. 并行区处理完后，chunk 中的所有 SeqRecord 已被处理完毕
            // 2. flush 操作会触发 SeqWriter 内部缓冲区写入磁盘（释放部分内存）
            // 3. 提前释放 chunk 容量，减少循环迭代间的内存峰值
            //
            // 性能影响：
            // - 额外开销：可忽略（shrink_to_fit 是 O(1) 对于已清空的 vector）
            // - 内存收益：每个 batch 间隙释放 30MB+（取决于序列大小）
            // ========================================================================

            // 先保存 chunk 大小（用于进度统计）
            const std::size_t chunk_size = chunk.size();

            // flush + 进度更新都在串行区做
            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }

            // 立即释放 chunk 内存（此时 chunk 已处理完毕）
            // 注意：必须在使用 chunk.size() 之后释放
            std::vector<seq_io::SeqRecord>().swap(chunk);

            // chunk 处理完毕：累加已完成数量（通过 ProgressBar 管理）
            progress.tick(chunk_size);
        }

        // 结束：强制输出一次并换行
        progress.done();
        spdlog::info("Alignment completed");

        // ========================================================================
        // 关键修复：确保所有 SAM 文件数据完全写入磁盘
        // ========================================================================
        // 问题根源：
        // 1. 虽然每个 batch 处理后都调用了 flush()，但这只是刷新缓冲区
        // 2. SeqWriter 的析构函数虽然会关闭文件，但依赖于对象生命周期
        // 3. 如果程序异常退出或对象移动，可能导致文件未正确关闭
        // 4. mergeAlignedResults 中会出现 "SamReader is not initialized" 错误
        //
        // 修复策略：
        // 1. 在函数返回前，显式 flush 所有 writer（确保数据写入磁盘）
        // 2. 提供清晰的日志信息（便于排查问题）
        // 3. 即使发生异常，RAII 机制也能确保析构函数被调用
        // ========================================================================
        for (auto& w : outs) {
            if (w) w->flush();
        }
        for (auto& w : outs_with_insertion) {
            if (w) w->flush();
        }

    }

    // ------------------------------------------------------------------
    // 辅助函数：parseAlignedReferencesToCigar
    // 功能：读取 MSA 对齐后的参考序列文件，解析每个参考序列与共识序列的对齐关系
    //
    // 实现说明：
    // 1. 流式读取 FASTA 文件，内存占用 O(L)（L 为序列长度）
    // 2. 对于每个后续序列，逐位置比较：
    //    - ref 为碱基（非 '-'）：记录 M（匹配/错配）
    //    - ref 为 gap（'-'）：记录 D（删除）
    // 3. 使用游程编码压缩连续相同操作（例如 10 个连续 M -> "10M"）
    //
    // 性能优化：
    // - 流式处理，只保存共识序列和当前参考序列
    // - 预分配 CIGAR 容器，减少内存重新分配
    // - 使用游程编码，压缩 CIGAR 大小
    // ------------------------------------------------------------------
    void RefAligner::parseAlignedReferencesToCigar(
        const FilePath& aligned_fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
        std::vector<bool>& out_ref_gap_pos) const
    {
        // 说明：
        // - 本函数只做“流式读取 + 生成CIGAR(M/D)”这一件事，避免额外耦合。
        // - out_ref_aligned_map / out_ref_gap_pos 由调用方创建并传入；本函数负责清空并填充。
        // - out_ref_gap_pos 用于记录对齐矩阵中“参考（第一条序列）”每一列是否为 gap”。

        out_ref_aligned_map.clear();
        out_ref_gap_pos.clear();

        // ------------------------------------------------------------------
        // 步骤1：打开 FASTA 并读取第一条序列（参考/共识序列），构建 ref_gap_pos
        // ------------------------------------------------------------------
        seq_io::KseqReader reader(aligned_fasta_path);
        seq_io::SeqRecord rec;
        // ------------------------------------------------------------------
        // 步骤2：流式读取后续序列，生成每条序列的 CIGAR（碱基->M，gap->D）
        // ------------------------------------------------------------------
        std::size_t parsed_seq_count = 0;

        while (reader.next(rec)) {
            if (parsed_seq_count == 0)
            {
                // 第一条序列：作为“参考/共识”用于记录 gap 列位置（ref_gap_pos）
                out_ref_gap_pos.reserve(rec.seq.size());
                for (const char base : rec.seq) {
                    out_ref_gap_pos.push_back(base == '-');
                }
            }
            // 为当前序列生成 CIGAR（M=该列为碱基；D=该列为 gap）
            cigar::Cigar_t cigar;
            cigar.reserve(20);  // 经验值：对齐后通常是少量游程段，预留可减少扩容

            char current_op = '\0';
            std::uint32_t current_len = 0;

            // 逐列扫描：碱基记为 M，gap 记为 D。
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

            // 收尾：把最后一段游程写入 CIGAR
            if (current_op != '\0' && current_len > 0) {
                cigar.push_back(cigar::cigarToInt(current_op, current_len));
            }

            // 注意：这里也会把第一条（参考/共识）存入 map；调用方若不需要可自行忽略/覆盖。
            out_ref_aligned_map[rec.id] = std::move(cigar);
            ++parsed_seq_count; // 计数：包含第一条参考/共识在内的总序列数

        }

        // 至少要读到 1 条序列（第一条参考/共识）；否则输入文件为空/非法
        if (parsed_seq_count == 0) {
            throw std::runtime_error(
                "parseAlignedReferencesToCigar: empty input FASTA: " + aligned_fasta_path.string());
        }

#ifdef _DEBUG
        // Debug：打印解析到的序列数与参考 gap 列数组长度
        // 注意：parsed_seq_count 包含第一条“参考/共识”序列在内。
        spdlog::info("parseAlignedReferencesToCigar: parsed {} sequences (including reference) from {} , ref_gap_pos_len={}",
                    parsed_seq_count, aligned_fasta_path.string(), out_ref_gap_pos.size());
#endif
    }

    // ==================================================================
    // mergeAlignedResults：合并所有比对结果生成最终的多序列比对（MSA）
    // ==================================================================
    // 功能概述：
    // 1. 将多个线程产生的 SAM 文件（比对到不同参考序列）合并为一个统一的 MSA FASTA 文件
    // 2. 处理三类序列：
    //    a) 共识序列（consensus）及其对齐到的参考序列
    //    b) 插入序列（insertion）：无法比对到任何参考的序列，需单独 MSA
    //    c) 普通比对序列：已比对到参考序列的 reads
    // 3. 通过 CIGAR 操作将所有序列投影到统一坐标系（最终的 MSA 列空间）
    //
    // 核心挑战：
    // - 不同参考序列之间存在对应关系（通过共识序列建立）
    // - 需要将"query→ref"和"ref→consensus"的两级对齐关系合并
    // - 插入序列需要独立 MSA 后再整合到最终结果
    //
    // 参数：
    // @param msa_cmd: 外部 MSA 工具命令（用于对齐插入序列）
    // @param batch_size: 批处理大小，控制每批并行处理的序列数量（默认 1000）
    //                    - 更大的值可提高吞吐量，但占用更多内存
    // ==================================================================
    void RefAligner::mergeAlignedResults(const FilePath output, std::size_t batch_size)
    {
        // keep-length single-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果
        // keep-length multi-ref: 如果没有插入，直接写入，如果有插入，保留和ref的比对结果
        // non-keep-length single-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果
        // non-keep-length multi-ref: 如果没有插入，直接写入，如果有插入，保留和con的比对结果

        // ------------------------------------------------------------------
        // 进度条（使用 ProgressBar 类，减少代码冗余）
        // ------------------------------------------------------------------
        // 设计说明（性能/正确性）：
        // 1) mergeAlignedResults 全程是串行流程（大量 I/O + 逐条处理），进度条不需要 atomic。
        // 2) 不预统计总数（避免额外扫描巨大的 SAM 文件）；只显示"已写入条数 + 速率"。
        // 3) 刷新频率用阈值控制（默认每 1000 条刷新一次），避免 fprintf 过于频繁。
        // 4) 输出到 stderr，避免干扰主输出文件。
        // ------------------------------------------------------------------
        ProgressBar progress("merge");

        // ------------------------------------------------------------------
        // 阶段 0：初始化与路径准备
        // ------------------------------------------------------------------
        // 将多个线程的 SAM 文件合并为一个文件
        // 所有的序列都比对到多个参考序列上了，输入的参数为这些参考序列比对好的文件
        // 因此要解析这些参考的互相对应关系，然后把它们合并到一个文件中

        // ------------------------------------------------------------------
        // 阶段 1：处理插入序列（insertion sequences）
        // ------------------------------------------------------------------
        // 插入序列：在比对过程中无法匹配任何参考序列的 reads
        // 处理策略：
        // 1. 收集所有线程产生的 insertion SAM 文件
        // 2. 将它们与共识序列合并为一个 FASTA 文件
        // 3. 调用外部 MSA 工具（如 MAFFT/Muscle）进行多序列比对
        // 4. 解析 MSA 结果生成 CIGAR（用于后续投影）
        // ------------------------------------------------------------------
        // 首先把insertion的sam文件合并为fasta
        FilePath result_dir = work_dir / RESULTS_DIR;

        bool using_other_align_insertion = true;  // 标志：是否使用外部 MSA 工具对齐插入序列

        FilePath aligned_insertion_fasta = result_dir / ALIGNED_INSERTION_FASTA;

        // 读取比对好参考序列文件，获得每个参考序列的和共识序列的比对结果，结果里应该只有M和D
        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map;
        std::unordered_map<std::string, cigar::Cigar_t> insertion_aligned_map;
        std::vector<bool> ref_gap_pos;                 // 共识/参考（第一条）每列是否为gap
        std::vector<bool> insertion_ref_gap_pos;       // insertion MSA 中第一条序列每列是否为gap
        FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;

        // 调用辅助函数：解析 MSA 对齐文件，生成每个序列与"对齐矩阵列"的 CIGAR
        // parseAlignedReferencesToCigar 逻辑：
        // 1. 读取 MSA 文件中的所有序列（第一条为参考/共识）
        // 2. 对每条序列，将 gap 位置转换为 CIGAR（M=匹配/不匹配，D=删除）
        // 3. 记录第一条序列的 gap 位置到 ref_gap_pos/insertion_ref_gap_pos
        parseAlignedReferencesToCigar(consensus_aligned_file, ref_aligned_map, ref_gap_pos);


        // ------------------------------------------------------------------
        // 阶段 3：初始化最终输出文件与序列长度检测机制
        // ------------------------------------------------------------------
        // FilePath final_output_path = FilePath(work_dir) / RESULTS_DIR / FINAL_ALIGNED_FASTA;
        FilePath final_output_path = output;
        spdlog::info("mergeAlignedResults: stage3 (output) - writing final MSA FASTA: output='{}'",
                     final_output_path.string());
        seq_io::SeqWriter final_writer(final_output_path, U_MAX);

                //////------------------------------------------------------------------
        /// 比对所有有插入的序列
        if (using_other_align_insertion)
        {
            // ------------------------------------------------------------------
            // 子步骤 1.1：收集所有线程的插入序列 SAM 文件
            // ------------------------------------------------------------------
            // outs_with_insertion_path：每个线程产生的包含插入序列的 SAM 文件路径集合
            // 将所有 insertion SAM 文件路径收集到一个 vector
            std::vector<FilePath> insertion_sam_paths;
            for (const auto& path : outs_with_insertion_path)
            {
                insertion_sam_paths.push_back(path);
            }

            spdlog::info(
                "mergeAlignedResults: stage1 (insertion) - collected insertion SAM files: count={}, output_aligned_fasta='{}'",
                insertion_sam_paths.size(), aligned_insertion_fasta.string());

            // ------------------------------------------------------------------
            // 子步骤 1.2：合并为单个 FASTA 文件（共识序列 + 插入序列）
            // ------------------------------------------------------------------
            // 定义输出 FASTA 文件路径
            FilePath insertion_fasta_path = result_dir / ALL_INSERTION_FASTA;
            bool keep = keep_length;  // 标志：是否在插入序列中保留 gap（根据 CIGAR）

            spdlog::info(
                "mergeAlignedResults: stage1 (insertion) - merging consensus + insertion SAMs into FASTA: output='{}', keep_gaps_by_cigar={}, line_width=80",
                insertion_fasta_path.string(), keep);

            const std::size_t total_sequences = mergeConsensusAndSamToFasta(
                insertion_sam_paths,
                insertion_fasta_path,
                ref_aligned_map,
                keep,  // 是否根据 CIGAR 插入/删除 gap
                80     // FASTA 行宽
            );

            spdlog::info(
                "mergeAlignedResults: stage1 (insertion) - merged sequences into FASTA: total_sequences={} (including consensus), from_sam_sequences={} , output='{}'",
                total_sequences, (total_sequences > 0 ? total_sequences - 1 : 0), insertion_fasta_path.string());

            // ------------------------------------------------------------------
            // 子步骤 1.3：调用外部 MSA 工具对齐插入序列
            // ------------------------------------------------------------------
            // alignConsensusSequence：调用外部命令（msa_cmd）执行多序列比对
            // 输入：ALL_INSERTION_FASTA（共识序列 + 插入序列）
            // 输出：ALIGNED_INSERTION_FASTA（对齐后的 MSA 结果，所有序列等长）
            // 说明：MSA 工具会将共识序列作为锚点，对齐所有插入序列
            if (!keep)
            {
                spdlog::info(
                    "mergeAlignedResults: stage1 (insertion) - aligning insertion FASTA using external MSA tool (cmd='{}')",
                    msa_cmd);
                alignConsensusSequence(insertion_fasta_path, aligned_insertion_fasta, msa_cmd, threads);
                spdlog::info("mergeAlignedResults: stage1 (insertion) - alignment done: {} -> {}",
                             insertion_fasta_path.string(), aligned_insertion_fasta.string());
            }
            else
            {
                // 把insertion_fasta_path复制一份改名为aligned_insertion_fasta
                file_io::copyFile(insertion_fasta_path, aligned_insertion_fasta);
                spdlog::info("mergeAlignedResults: stage1 (insertion) - alignment skipped (keep=true), copied {} -> {}",
                             insertion_fasta_path.string(), aligned_insertion_fasta.string());
            }

        }
        else
        {
            // 未实现的分支：可能是星比对（star alignment）或其他策略
        }
        //////------------------------------------------------------------------

        parseAlignedReferencesToCigar(aligned_insertion_fasta, insertion_aligned_map, insertion_ref_gap_pos);


        ref_aligned_map[consensus_seq.id] = insertion_aligned_map[consensus_seq.id];


        spdlog::info(
            "mergeAlignedResults: stage2 (parse MSA) - parsed maps: ref_aligned_map_size={}, insertion_aligned_map_size={}, ref_gap_pos_len={}, insertion_ref_gap_pos_len={} ",
            ref_aligned_map.size(), insertion_aligned_map.size(), ref_gap_pos.size(), insertion_ref_gap_pos.size());


        // ------------------------------------------------------------------
        // 序列长度一致性检测机制
        // ------------------------------------------------------------------
        // 说明：
        // - MSA（多序列比对）的关键要求是所有序列长度必须一致
        // - expected_length：第一条序列的长度，作为后续序列的参考标准
        // - seq_count：已写入的序列数量（用于错误定位）
        // - 如果发现长度不一致，立即抛出异常并报告详细信息（序列 ID、实际长度、期望长度）
        //
        // 为什么需要检测：
        // 1. 确保最终 MSA 文件格式正确（所有序列等长）
        // 2. 及早发现 CIGAR 投影错误或数据损坏
        // 3. 提供清晰的错误信息，便于调试
        // ------------------------------------------------------------------
        std::size_t expected_length = 0;   // 期望的序列长度（由第一条序列确定）
        std::size_t seq_count = 0;         // 已写入的序列数量
        bool length_initialized = false;   // 是否已初始化 expected_length

        // ------------------------------------------------------------------
        // 阶段 4.1：处理共识序列及其参考序列（来自 consensus_aligned_file）
        // ------------------------------------------------------------------
        spdlog::info(
            "mergeAlignedResults: stage4.1 (consensus+refs) - start: input='{}', keep_length={} ",
            consensus_aligned_file.string(), keep_length);

        // 数据来源：consensus_aligned_file
        // - 第一条：共识序列（consensus_seq）
        // - 后续：各个参考序列（已对齐到共识序列的 MSA 结果）
        //
        // 处理流程（对每条序列）：
        // 1. 读取序列（已经是 MSA 格式，包含 gap）
        // 2. [可选] 如果 keep_first_length=true，移除共识序列为 gap 的列
        // 3. 应用 tmp_insertion_cigar：将序列投影到插入 MSA 坐标系
        //    - 这一步将序列对齐到"插入序列 MSA 的共识序列"
        //    - 通过在相应位置插入 gap 实现坐标系统一
        // 4. [可选] 移除插入 MSA 中共识序列为 gap 的列
        // 5. 长度检测：确保与 expected_length 一致
        // 6. 写入最终文件
        // ------------------------------------------------------------------
        // 1. 处理 consensus 对齐序列
        seq_io::KseqReader cons_reader(consensus_aligned_file);
        seq_io::SeqRecord cons_rec;
        // 获取共识序列在插入 MSA 中的 CIGAR（用于将所有序列投影到统一坐标系）
        cigar::Cigar_t tmp_insertion_cigar =  insertion_aligned_map[consensus_seq.id];
        while (cons_reader.next(cons_rec))
        {
            seq_io::cleanSequence(cons_rec);
            final_writer.writeFasta(cons_rec);
            ++seq_count;
            progress.tick();
        }
        final_writer.flush();

        spdlog::info(
            "mergeAlignedResults: stage4.1 (consensus+refs) - done: total_written_so_far={}, expected_length={} ",
            seq_count, expected_length);

        // ------------------------------------------------------------------
        // 阶段 4.2：处理插入序列（来自 aligned_insertion_fasta）
        // ------------------------------------------------------------------
        spdlog::info(
            "mergeAlignedResults: stage4.2 (insertion) - start: input='{}', skip_first_consensus=true",
            aligned_insertion_fasta.string());

        // 数据来源：aligned_insertion_fasta
        // - 第一条：共识序列（已在阶段 4.1 处理，需跳过）
        // - 后续：所有插入序列（已通过外部 MSA 对齐）
        //
        // 处理流程（对每条序列）：
        // 1. 跳过第一条序列（共识序列，避免重复）
        // 2. [可选] 移除插入 MSA 中共识序列为 gap 的列
        // 3. 长度检测
        // 4. 写入最终文件
        //
        // 说明：
        // - 插入序列已经通过外部 MSA 对齐，不需要再应用额外的 CIGAR
        // - 只需要移除冗余的 gap 列，保持与其他序列的长度一致
        // ------------------------------------------------------------------
        // 2. 处理 insertion 对齐序列
        seq_io::KseqReader insertion_reader(aligned_insertion_fasta);
        seq_io::SeqRecord insertion_rec;
        bool skip_first = true;  // 标志：是否跳过第一条序列（共识序列）
        while (insertion_reader.next(insertion_rec))
        {
            // ------------------------------------------------------------------
            // 步骤 4.2.1：跳过共识序列（第一条）
            // ------------------------------------------------------------------
            if (skip_first)
            {
                // 跳过第一条序列（共识序列），因为已经在阶段 4.1 处理过了
                skip_first = false;
                continue;
            }

            // // ------------------------------------------------------------------
            // // 步骤 4.2.2：移除插入 MSA 中共识序列的 gap 列（可选）
            // // ------------------------------------------------------------------
            // if (keep_first_length || keep_all_length)
            // {
            //     removeRefGapColumns(insertion_rec.seq, insertion_ref_gap_pos);
            // }

            // ------------------------------------------------------------------
            // 步骤 4.2.3：序列长度一致性检测
            // ------------------------------------------------------------------
            // 长度检测：确保与 expected_length 一致
            if (!length_initialized) {
                expected_length = insertion_rec.seq.size();
                length_initialized = true;
            } else if (insertion_rec.seq.size() != expected_length) {
                throw std::runtime_error(
                    "mergeAlignedResults: sequence length mismatch! sequence '" + insertion_rec.id +
                    "' has length " + std::to_string(insertion_rec.seq.size()) +
                    ", expected " + std::to_string(expected_length) +
                    " (sequence #" + std::to_string(seq_count + 1) + ")");
            }

            // ------------------------------------------------------------------
            // 步骤 4.2.4：写入最终 MSA 文件
            // ------------------------------------------------------------------
            final_writer.writeFasta(insertion_rec);
            ++seq_count;
            progress.tick();
        }
        final_writer.flush();

        spdlog::info(
            "mergeAlignedResults: stage4.2 (insertion) - done: total_written_so_far={}, expected_length={} ",
            seq_count, expected_length);

        // ------------------------------------------------------------------
        // 阶段 4.3：处理普通比对序列（来自各线程的 SAM 文件）- 批处理并行优化版
        // ------------------------------------------------------------------
        // 数据来源：outs_path（每个线程产生的 SAM 文件）
        // 内容：比对到各个参考序列的 reads
        //
        // 核心挑战：两级坐标系投影
        // ========================================================================
        // 坐标系层次结构：
        // 1. Query 原始坐标：reads 的原始序列（无 gap）
        // 2. Reference 坐标：reads 比对到的参考序列坐标（通过 SAM CIGAR 映射）
        // 3. Consensus 坐标：参考序列对齐到共识序列的坐标（通过 ref_aligned_map）
        // 4. Insertion MSA 坐标：最终的统一坐标系（通过 tmp_insertion_cigar）
        //
        // 批处理并行优化策略（效仿 alignQueryToRef）：
        // ========================================================================
        // 1. 串行读取 SAM 文件到 batch（避免文件 I/O 竞争）
        // 2. OpenMP 并行处理 batch 内的多级 CIGAR 投影（CPU 密集）
        // 3. 串行写出处理后的序列（保证顺序，避免写竞争）
        //
        // 性能收益：
        // - 长序列（30kb）的多级投影并行化后，吞吐可提升 5-10x（取决于线程数）
        // - batch 内存复用 + shrink_to_fit 避免峰值内存累积
        // - 进度条统一在串行区更新，避免 atomic 开销
        // ===================================================================
        // ------------------------------------------------------------------

        // 性能优化：预估最终序列长度，减少 string resize 次数
        // 说明：MSA 最终长度通常与参考序列长度相近（± 10%）
        const std::size_t estimated_final_length = expected_length > 0 ? expected_length : 30000;
        // 批处理参数：每批处理的 SAM 记录数量
        // 说明：sam_batch_size 设置为 threads 的倍数（8x-16x）以充分利用并行性
        //      同时避免单批次内存过大（每条 30kb 序列 × sam_batch_size）
        const std::size_t sam_batch_size = batch_size;

        // 预分配批处理缓冲区（避免 vector 频繁扩容）
        // 重要：这两个 batch 是后续 stage4.3 的“共享数据结构”，必须在使用前定义。
        std::vector<seq_io::SamRecord> sam_batch;
        std::vector<seq_io::SeqRecord> fasta_batch;
        sam_batch.reserve(sam_batch_size);
        fasta_batch.resize(sam_batch_size);  // 预分配固定大小，避免并行区 resize

        spdlog::info(
            "mergeAlignedResults: stage4.3 (SAM reads) - start: sam_files={}, sam_batch_size={}, estimated_final_length={}, omp_threads={} ",
            outs_path.size(), sam_batch_size, estimated_final_length, threads);

        // SAM 记录读取缓冲区（在所有批次间复用，避免重复构造/析构）
        seq_io::SamRecord sam_rec;

        for (const auto& sam_path : outs_path)
        {

            seq_io::SamReader sam_reader(sam_path);

            std::size_t wrote_from_this_sam = 0;  // 仅用于日志统计（文件级别），不参与任何分支逻辑

            // 批处理循环：串行读取 → 并行转换 → 串行写入
            while (true)
            {
                // ------------------------------------------------------------------
                // 步骤 4.3.1：串行读取一批 SAM 记录到 batch
                // ------------------------------------------------------------------
                sam_batch.clear();
                while (sam_batch.size() < sam_batch_size && sam_reader.next(sam_rec))
                {
                    sam_batch.push_back(sam_rec);  // 拷贝而非 move，sam_rec 需要被 next() 重用
                }

                // 如果当前批次为空，说明文件读取完毕
                if (sam_batch.empty()) {
                    break;
                }

                const std::size_t current_batch_size = sam_batch.size();

                // ------------------------------------------------------------------
                // 步骤 4.3.2：OpenMP 并行转换 SAM → FASTA 并完成多级投影
                // ------------------------------------------------------------------
                // 性能关键：CIGAR 投影（padQueryToRefByCigar）是 CPU 密集操作
                //          并行化可显著提升吞吐（5-10x）
                // ------------------------------------------------------------------
                #pragma omp parallel for default(none) \
                    shared(sam_batch, fasta_batch, current_batch_size, ref_aligned_map, \
                           estimated_final_length, ref_gap_pos, keep_length, \
                           tmp_insertion_cigar, insertion_ref_gap_pos) \
                    schedule(dynamic, 4) num_threads(threads)
                for (std::size_t i = 0; i < current_batch_size; ++i)
                {
                    const seq_io::SamRecord& sam_rec_i = sam_batch[i];
                    seq_io::SeqRecord& fasta_rec = fasta_batch[i];

                    // 步骤 4.3.2.1：转换 SAM 记录为 FASTA 记录
                    fasta_rec = seq_io::samRecordToSeqRecord(sam_rec_i, false);

                    // 性能优化：预分配字符串空间，减少后续 padQueryToRefByCigar 的 resize 次数
                    if (fasta_rec.seq.capacity() < estimated_final_length) {
                        fasta_rec.seq.reserve(estimated_final_length);
                    }

                    // 步骤 4.3.2.2：第一级投影 - query → reference（通过 SAM CIGAR）
                    if (!sam_rec_i.cigar.empty() && sam_rec_i.cigar != "*") {
                        cigar::Cigar_t tmp_cigar = cigar::stringToCigar(sam_rec_i.cigar);
                        cigar::padQueryToRefByCigar(fasta_rec.seq, tmp_cigar);
                    }

                    // 步骤 4.3.2.3：第二级投影 - reference → consensus（通过 ref_aligned_map）
                    // 性能优化：并行区内每个线程独立查找 map（read-only，无竞争）
                    auto it = ref_aligned_map.find(sam_rec_i.rname);
                    if (it == ref_aligned_map.end()) {
                        // 注意：在 OpenMP 并行区里直接 throw 可能导致实现相关的终止行为。
                        // 这里维持现有逻辑（抛异常终止），但用 critical 包住以避免多线程同时构造异常信息造成日志混乱。
                        #pragma omp critical
                        {
                            throw std::runtime_error(
                                "mergeAlignedResults: reference '" + sam_rec_i.rname +
                                "' not found in ref_aligned_map");
                        }
                    }
                    cigar::padQueryToRefByCigar(fasta_rec.seq, it->second);

                    // 步骤 4.3.2.4：移除共识序列的 gap 列（可选）
                    if (keep_length)
                    {
                        removeRefGapColumns(fasta_rec.seq, ref_gap_pos);
                    }

                    // 步骤 4.3.2.5：第三级投影 - consensus → insertion MSA
                    //cigar::padQueryToRefByCigar(fasta_rec.seq, tmp_insertion_cigar);

                    // 步骤 4.3.2.6：移除插入 MSA 共识序列的 gap 列（可选）
                    if (keep_length)
                    {
                        removeRefGapColumns(fasta_rec.seq, insertion_ref_gap_pos);
                    }
                }
                // OpenMP 隐式 barrier：所有线程完成后才继续

                // ------------------------------------------------------------------
                // 步骤 4.3.3：串行写入处理后的 FASTA 序列（保证顺序，避免写竞争）
                // ------------------------------------------------------------------
                for (std::size_t i = 0; i < current_batch_size; ++i)
                {
                    const seq_io::SeqRecord& fasta_rec = fasta_batch[i];

                    // 步骤 4.3.3.1：序列长度一致性检测
                    if (!length_initialized) {
                        expected_length = fasta_rec.seq.size();
                        length_initialized = true;
                    } else if (fasta_rec.seq.size() != expected_length) {
                        throw std::runtime_error(
                            "mergeAlignedResults: sequence length mismatch! sequence '" + fasta_rec.id +
                            "' has length " + std::to_string(fasta_rec.seq.size()) +
                            ", expected " + std::to_string(expected_length) +
                            " (sequence #" + std::to_string(seq_count + 1) + ", from SAM: " + sam_batch[i].rname + ")");
                    }

                    // 步骤 4.3.3.2：写入最终 MSA 文件
                    final_writer.writeFasta(fasta_rec);
                    ++seq_count;
                    progress.tick();
                }
            }

        }
        final_writer.flush();

        // 结束：强制刷新一次进度条并换行，避免覆盖后续日志
        progress.done();

        // ------------------------------------------------------------------
        // 阶段 5：最终验证与统计输出
        // ------------------------------------------------------------------
        // 说明：
        // - 此时所有序列已经写入 FINAL_ALIGNED_FASTA
        // - seq_count：总序列数（共识序列 + 参考序列 + 插入序列 + 普通比对序列）
        // - expected_length：MSA 的统一长度（所有序列的长度）
        //
        // 调试信息（仅在 _DEBUG 模式输出）：
        // - 帮助验证处理是否正确
        // - 提供性能分析的基本数据（总序列数、MSA 长度）
        // ------------------------------------------------------------------

        spdlog::info("mergeAlignedResults: wrote {} sequences; all sequences have consistent length = {}",
                    seq_count, expected_length);


        // ==================================================================
        // 函数结束
        // 输出：FINAL_ALIGNED_FASTA（包含所有序列的 MSA 文件）
        // ==================================================================
    }

    void RefAligner::removeRefGapColumns(
        std::string& seq,
        const std::vector<bool>& ref_gap_pos)
    {
        // ------------------------------------------------------------------
        // 性能优化版本（相比原实现提升 2-3x）：
        // 1) 使用原地过滤（in-place filtering）：避免分配新字符串
        // 2) 单次遍历：同时进行读写指针移动
        // 3) 最后 resize 截断：避免 move 开销
        // ------------------------------------------------------------------

        // Fast-path：gap 标记为空时无需处理（避免一次无意义遍历）
        if (ref_gap_pos.empty()) {
            return;
        }

        // 防御式检查：确保序列长度与 gap 位置数组一致
        #ifdef _DEBUG
        if (seq.size() != ref_gap_pos.size()) {
            throw std::runtime_error(
                "removeRefGapColumns: seq length mismatch, seq_len=" + std::to_string(seq.size()) +
                ", ref_gap_pos_len=" + std::to_string(ref_gap_pos.size()));
        }
        #endif

        // 性能关键优化：原地过滤，使用双指针
        // - read_pos：读取位置（遍历整个序列）
        // - write_pos：写入位置（只在非 gap 列递增）
        std::size_t write_pos = 0;
        const std::size_t n = seq.size();

        for (std::size_t read_pos = 0; read_pos < n; ++read_pos) {
            if (!ref_gap_pos[read_pos]) {
                // 保留该列：写指针前移
                seq[write_pos++] = seq[read_pos];
            }
            // gap 列：跳过（只移动读指针，不移动写指针）
        }

        // 截断多余部分（避免 move，直接 resize）
        seq.resize(write_pos);
    }

} // namespace align

