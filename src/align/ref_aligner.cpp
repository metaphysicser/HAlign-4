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

    namespace {
        std::string profileToGappedSequence(const ProfileMatrix& profile)
        {
            if (profile.len < 0 || profile.dim < 5) {
                throw std::runtime_error("profileToGappedSequence: invalid profile shape");
            }

            const std::size_t len = static_cast<std::size_t>(profile.len);
            const std::size_t dim = static_cast<std::size_t>(profile.dim);
            if (profile.prof.size() < len * dim) {
                throw std::runtime_error("profileToGappedSequence: profile data is truncated");
            }

            static constexpr char bases[5] = {'A', 'C', 'G', 'T', 'N'};
            std::string seq;
            seq.reserve(len);

            for (std::size_t col = 0; col < len; ++col) {
                const std::size_t off = col * dim;
                std::uint64_t total = 0;
                std::uint32_t best_count = 0;
                std::size_t best_idx = 4;

                for (std::size_t idx = 0; idx < 5; ++idx) {
                    const std::uint32_t count = profile.prof[off + idx];
                    total += count;
                    if (count > best_count) {
                        best_count = count;
                        best_idx = idx;
                    }
                }

                seq.push_back(total == 0 ? '-' : bases[best_idx]);
            }

            return seq;
        }
    } // namespace

    // 读取参考序列，计算索引，生成共识序列
        RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                                                     int kmer_size, int window_size,
                                                                                                         int sketch_size, int sketch_kmer_size, bool noncanonical,
                                                     int threads, std::string msa_cmd,
                                                     bool keep_length,
                                                     bool enable_wfa,
                                                     const FilePath& ref_aligned_path,
                                                     std::array<int8_t, 25> score_matrix_in,
                                                     int gap_open,
                                                     int gap_extend)
        : work_dir(work_dir),
          kmer_size(kmer_size),
          window_size(window_size),
          sketch_size(sketch_size),
		  sketch_kmer_size(sketch_kmer_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_cmd(std::move(msa_cmd)),
          keep_length(keep_length),
          enable_wfa(enable_wfa),
          score_matrix(score_matrix_in),
          gap_open(gap_open),
          gap_extend(gap_extend)
    {
	        // 加载参考序列并构建 sketch/minimizer 索引
	        seq_io::KseqReader reader(ref_fasta_path);
	        seq_io::SeqRecord rec;
        	spdlog::info("Loading reference sequences from {}", ref_fasta_path.string());

	        while (reader.next(rec)) {
	            // 关键改动：sketch 使用独立的 sketch_kmer_size，
	            // minimizer 仍使用 kmer_size，保持锚点密度/行为不变。
	            auto sketch = mash::sketchFromSequence(rec.seq, sketch_kmer_size, sketch_size,
	                                                   noncanonical, random_seed);
	        	ref_sequences.insert({rec.id, rec});
	            ref_sketch.insert({rec.id, std::move(sketch)});
	        }
        	spdlog::info("Loaded {} reference sequences", ref_sequences.size());

	        // 设置共识序列生成的文件路径
	        const bool has_prealigned_ref = !ref_aligned_path.empty();
	        const FilePath consensus_unaligned_file = has_prealigned_ref ? ref_aligned_path : ref_fasta_path;
	        const FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
	        const FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
	        const FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

	        // 执行 MSA 并生成共识序列
	        // 说明：
	        // - 默认路径：对 -r/--ref 提供的原始 FASTA 重新做一次 MSA，再从对齐结果生成共识；
	        // - 若用户同时提供 --ref-align，则说明这个参考集合已经有现成的 MSA，
	        //   这里直接复用该 MSA，避免重复对齐，从而节省时间并保持参考坐标不变。
	        constexpr std::size_t consensus_batch_size = 4096;
	        if (has_prealigned_ref) {
	            file_io::copyFile(consensus_unaligned_file, consensus_aligned_file);
	        } else {
	            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, this->msa_cmd, threads);
	        }

        	seq_io::KseqReader reader2(consensus_aligned_file);
        	seq_io::SeqRecord rec2;
        	while (reader2.next(rec2))
        	{
        		ref_profile.insert({rec2.id, ProfileMatrix::fromAlignedSequences({rec2.seq})});
        	}
	        consensus::ConsensusResult consensus_result = consensus::generateConsensusResult(
	            consensus_aligned_file, consensus_file, consensus_json_file,
	            0, threads, consensus_batch_size);

	        consensus_gap_seq.id = "consensus";
	        consensus_gap_seq.seq = std::move(consensus_result.gap_seq);
	        consensus_profile = ProfileMatrix::fromConsensusCounts(consensus_result.counts);

        	// 以下变量为seq2seq准备
        	consensus_seq.id = "consensus";
        	consensus_seq.seq = std::move(consensus_result.seq);
	        // 预计算共识序列的 sketch 和 minimizer，避免重复计算
	        consensus_sketch = mash::sketchFromSequence(
	            consensus_seq.seq,
	            static_cast<std::size_t>(sketch_kmer_size),
	            static_cast<std::size_t>(sketch_size),
	            noncanonical,
	            random_seed);

	        consensus_minimizer = minimizer::extractMinimizer(
	            consensus_seq.seq, kmer_size, window_size, noncanonical);
	        consensus_gap_minimizer = minimizer::extractMinimizer(
	            consensus_gap_seq.seq, kmer_size, window_size, noncanonical);
	    }

    // Options 配置委托构造
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,
            ref_fasta_path,
            opt.kmer_size,
            opt.kmer_window,
            opt.sketch_size,
            opt.sketch_kmer_size,
            true,
            opt.threads,
            opt.msa_cmd,
            opt.keep_length,
            opt.wfa,
            FilePath(opt.ref_align_path),
            opt.score_matrix,
            opt.gap_open,
            opt.gap_extend)
    {
    }

    AlignConfig RefAligner::makeAlignConfig() const
    {
        AlignConfig cfg;
        cfg.mat = score_matrix.data();
        cfg.gap_open = gap_open;
        cfg.gap_extend = gap_extend;
        return cfg;
    }

    // 全局比对：生成 minimizer 锚点，执行比对
    cigar::Cigar_t RefAligner::Seq2SeqWithAnchor(const std::string& ref,
                                           const std::string& query,
                                           double similarity,
                                           const SeedHits* ref_minimizer,
                                           const SeedHits* query_minimizer) const
    {
        const SeedHits* ref_mz_ptr = ref_minimizer;
        const SeedHits* qry_mz_ptr = query_minimizer;

        SeedHits ref_mz_tmp;
        SeedHits qry_mz_tmp;

        // 若 minimizer 为空，现场计算
        if (ref_mz_ptr == nullptr || ref_mz_ptr->empty()) {
            ref_mz_tmp = minimizer::extractMinimizer(ref, kmer_size, window_size, noncanonical);
            ref_mz_ptr = &ref_mz_tmp;
        }
        if (qry_mz_ptr == nullptr || qry_mz_ptr->empty()) {
            qry_mz_tmp = minimizer::extractMinimizer(query, kmer_size, window_size, noncanonical);
            qry_mz_ptr = &qry_mz_tmp;
        }

        const anchor::Anchors anchors = minimizer::collect_anchors(*ref_mz_ptr, *qry_mz_ptr);
        cigar::Cigar_t result = globalAlignSeq2Seq(ref, query, anchors, makeAlignConfig());

#ifdef _DEBUG
        const std::size_t cigar_ref_len = cigar::getRefLength(result);
        const std::size_t cigar_qry_len = cigar::getQueryLength(result);
        if (cigar_ref_len != ref.size() || cigar_qry_len != query.size()) {
            spdlog::debug("globalAlign: CIGAR length mismatch! ref:{} vs {}, query:{} vs {}",
                         ref.size(), cigar_ref_len, query.size(), cigar_qry_len);
        }
#endif
        return result;
    }

    cigar::Cigar_t RefAligner::Seq2ProfileWithAnchor(const ProfileMatrix& ref,
                                        const std::string& ref_string,
                                       const std::string& query,
                                       double similarity,
                                       int thread,
                                       const SeedHits* ref_minimizer,
                                       const SeedHits* query_minimizer) const
    {
        const SeedHits* ref_mz_ptr = ref_minimizer;
        const SeedHits* qry_mz_ptr = query_minimizer;

        SeedHits ref_mz_tmp;
        SeedHits qry_mz_tmp;

        // 若 minimizer 为空，现场计算
        if (ref_mz_ptr == nullptr || ref_mz_ptr->empty()) {
            ref_mz_tmp = minimizer::extractMinimizer(ref_string, kmer_size, window_size, noncanonical);
            ref_mz_ptr = &ref_mz_tmp;
        }
        if (qry_mz_ptr == nullptr || qry_mz_ptr->empty()) {
            qry_mz_tmp = minimizer::extractMinimizer(query, kmer_size, window_size, noncanonical);
            qry_mz_ptr = &qry_mz_tmp;
        }

        const anchor::Anchors anchors = minimizer::collect_anchors(*ref_mz_ptr, *qry_mz_ptr);
            cigar::Cigar_t result;
        if (thread > 0){
            result = globalAlignSeq2ProfileParallel(ref, ref_string, query, anchors, thread, makeAlignConfig());
        } else {
            result = globalAlignSeq2Profile(ref, ref_string, query, anchors, makeAlignConfig());
        }
#ifdef _DEBUG
        const std::size_t cigar_ref_len = cigar::getRefLength(result);
        const std::size_t cigar_qry_len = cigar::getQueryLength(result);
        if (cigar_ref_len != ref_string.size() || cigar_qry_len != query.size()) {
            spdlog::debug("globalAlign: CIGAR length mismatch! ref:{} vs {}, query:{} vs {}",
                         ref_string.size(), cigar_ref_len, query.size(), cigar_qry_len);
        }
#endif
        return result;

    }

    // 写入 SAM 记录
    void RefAligner::writeSamRecord(const seq_io::SeqRecord& q,
                                    const cigar::Cigar_t& cigar,
                                    std::string_view ref_name,
                                    seq_io::SeqWriter& out) const
    {
        const std::string cigar_str = cigar::cigarToString(cigar);
        const auto sam_rec = seq_io::makeSamRecord(q, ref_name, cigar_str, 1, 60, 0);
        out.writeSam(sam_rec);
    }

    // 单条 query 比对
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
    {
        // 计算 query 的 sketch 和 minimizer
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(sketch_kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            q.seq, kmer_size, window_size, noncanonical);

        const double consensus_similarity = mash::jaccard(qsk, consensus_sketch);
        cigar::Cigar_t consensus_cigar = Seq2SeqWithAnchor(
            consensus_seq.seq, q.seq, consensus_similarity,
            &consensus_minimizer, &query_minimizer);

        if (cigar::hasInsertion(consensus_cigar)) {
            writeSamRecord(q, consensus_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, consensus_cigar, consensus_seq.id, out);
        }
    }

    void RefAligner::alignOneQueryToProfile(const seq_io::SeqRecord& q,
                       seq_io::SeqWriter& out,
                       seq_io::SeqWriter& out_insertion,
                       cigar::Cigar_t& out_cigar,
                       std::string& out_ref_id,
                       int thread) const
    {
        // 初始化输出占位，便于调用方在调试时识别“未写入”状态
        out_cigar.clear();
        out_ref_id.clear();

        // 计算 query 的 sketch 和 minimizer
        const mash::Sketch qsk = mash::sketchFromSequence(
            q.seq,
            static_cast<std::size_t>(sketch_kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            q.seq, kmer_size, window_size, noncanonical);

        const ProfileMatrix* best_ref_profile = &consensus_profile;
        std::string best_ref_string = consensus_gap_seq.seq;
        std::string best_ref_id;
        const SeedHits* best_ref_minimizer = nullptr;
        SeedHits best_ref_minimizer_storage;
        double best_jaccard = -1.0;

        if (ref_sequences.size() > 1) {
            for (const auto& [ref_id, sketch] : ref_sketch) {
                const double j = mash::jaccard(qsk, sketch);
                if (best_ref_id.empty() || j > best_jaccard ||
                    (j == best_jaccard && ref_id < best_ref_id)) {
                    best_jaccard = j;
                    best_ref_id = ref_id;
                }
            }

            auto profile_it = ref_profile.find(best_ref_id);
            if (profile_it == ref_profile.end()) {
                throw std::runtime_error("Reference profile '" + best_ref_id + "' not found");
            }

            best_ref_profile = &profile_it->second;
            best_ref_string = profileToGappedSequence(*best_ref_profile);
            best_ref_minimizer_storage = minimizer::extractMinimizer(
                best_ref_string, kmer_size, window_size, noncanonical);
            best_ref_minimizer = &best_ref_minimizer_storage;
        } else {
            best_ref_minimizer = &consensus_gap_minimizer;
            best_jaccard = mash::jaccard(qsk, consensus_sketch);
        }

        // 执行全局比对
        cigar::Cigar_t initial_cigar = Seq2ProfileWithAnchor(
            *best_ref_profile, best_ref_string, q.seq, best_jaccard, thread,
            best_ref_minimizer, &query_minimizer);

        // 单参考序列：直接使用初始比对结果
        if (ref_sequences.size() == 1) {
            // 记录“最终输出”的 CIGAR 和参考 id，便于批处理阶段复用/统计
            out_cigar = initial_cigar;
            out_ref_id.clear();

            if (cigar::hasInsertion(initial_cigar)) {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out_insertion);
            } else {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out);
            }
            return;
        }

        // 多参考序列 + keep_length：保留与最佳参考的比对结果
        if (keep_length) {
            out_cigar = initial_cigar;
            out_ref_id = best_ref_id;

            if (cigar::hasInsertion(initial_cigar)) {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out_insertion);
            } else {
                writeSamRecord(q, initial_cigar, consensus_seq.id, out);
            }
            return;
        }

        // 多参考序列 + 非 keep_length：有插入时与共识序列二次比对
        const double consensus_similarity = mash::jaccard(qsk, consensus_sketch);
        cigar::Cigar_t recheck_cigar = Seq2ProfileWithAnchor(
            consensus_profile, consensus_gap_seq.seq, q.seq, consensus_similarity, thread,
            &consensus_gap_minimizer, &query_minimizer);

        out_cigar = recheck_cigar;
        out_ref_id.clear(); // 使用共识作为最终参考，避免误解为某条原始参考序列

        if (cigar::hasInsertion(recheck_cigar)) {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(q, recheck_cigar, consensus_seq.id, out);
        }
    }

    bool RefAligner::applyCigarToProfile(
        const std::string& query_seq,
        const cigar::Cigar_t& cigar,
        ProfileMatrix& target_profile)
    {
        // 关键约束：这里只做“计数累加”，不改变 profile 的列数/维度，避免影响现有比对逻辑。
        const std::size_t profile_len = static_cast<std::size_t>(target_profile.len);
        const std::size_t profile_dim = static_cast<std::size_t>(target_profile.dim);

        if (profile_dim == 0 || target_profile.prof.size() != profile_len * profile_dim) {
            return false;
        }

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        for (const cigar::CigarUnit unit : cigar) {
            char op = '\0';
            std::uint32_t len = 0;
            cigar::intToCigar(unit, op, len);

            switch (op) {
                case 'M':
                case '=':
                case 'X': {
                    // M/=/X 同时消耗参考和 query：把 query 当前碱基累加到对应参考列。
                    for (std::uint32_t k = 0; k < len; ++k) {
                        if (ref_pos >= profile_len || qry_pos >= query_seq.size()) {
                            return false;
                        }
                        const std::size_t base_idx = static_cast<std::size_t>(
                            align::ScoreChar2Idx[static_cast<unsigned char>(query_seq[qry_pos])]);
                        ++target_profile.prof[ref_pos * profile_dim + base_idx];
                        ++ref_pos;
                        ++qry_pos;
                    }
                    break;
                }
                case 'D':
                case 'N': {
                    // D/N 只消耗参考：该列没有 query 碱基贡献，仅推进参考坐标。
                    ref_pos += static_cast<std::size_t>(len);
                    if (ref_pos > profile_len) {
                        return false;
                    }
                    break;
                }
                case 'I':
                case 'S': {
                    // I/S 只消耗 query：不对应参考列，不能写入 profile，直接推进 query 坐标。
                    qry_pos += static_cast<std::size_t>(len);
                    if (qry_pos > query_seq.size()) {
                        return false;
                    }
                    break;
                }
                case 'H':
                case 'P': {
                    // H/P 不消耗 query 序列字符串内容，也不消耗参考列，保持坐标不变。
                    break;
                }
                default:
                    return false;
            }
        }

        // 与目标 profile 对齐时，参考消耗长度必须精确覆盖 profile 全长，避免越界和错列更新。
        return ref_pos == profile_len;
    }

    void RefAligner::updateProfilesFromChunk(
        const std::vector<seq_io::SeqRecord>& chunk,
        const std::vector<cigar::Cigar_t>& cigar_chunk,
        const std::vector<std::string>& ref_id_chunk)
    {
        if (chunk.size() != cigar_chunk.size() || chunk.size() != ref_id_chunk.size()) {
            throw std::runtime_error("updateProfilesFromChunk: chunk/cigar/ref_id size mismatch");
        }

        // 串行更新共享 profile：避免在并行区对同一 profile 加锁，减少锁竞争与缓存抖动。
        for (std::size_t i = 0; i < chunk.size(); ++i) {
            if (cigar_chunk[i].empty()) {
                continue;
            }

            ProfileMatrix* target_profile = nullptr;
            if (ref_id_chunk[i].empty()) {
                // 空 id 表示最终参考为共识序列，更新 consensus_profile。
                target_profile = &consensus_profile;
            } else {
                auto profile_it = ref_profile.find(ref_id_chunk[i]);
                if (profile_it != ref_profile.end()) {
                    target_profile = &profile_it->second;
                } else {
#ifdef _DEBUG
                    spdlog::debug("updateProfilesFromChunk: skip invalid ref_id={} at i={}",
                                  ref_id_chunk[i], i);
#endif
                    continue;
                }
            }

            if (applyCigarToProfile(chunk[i].seq, cigar_chunk[i], *target_profile)) {
                // depth 表示累计纳入 profile 的序列数，成功更新一条后再增加。
                ++target_profile->depth;
            } else {
#ifdef _DEBUG
                spdlog::debug("updateProfilesFromChunk: skip invalid cigar for query={} at i={}",
                              chunk[i].id, i);
#endif
            }
        }
    }

    // 批量比对 query 序列 - 并行处理，每线程独立输出
    void RefAligner::alignSeq2Seq(const FilePath& qry_fasta_path, std::size_t batch_size)
    {
        profile_alignment_mode = false;

        // 参数检查和初始化
        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: reference sequence is empty");
        }

        constexpr std::size_t default_batch_size = 2560;
        if (batch_size == 0) {
            batch_size = default_batch_size;
        }

        // 设置线程数
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        const FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        spdlog::info("Starting alignment: {} threads, batch size {}", nthreads, batch_size);

        // 为每个线程创建独立的输出文件
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

        // 流式读取 + 批处理并行
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        chunk.reserve(batch_size);

        ProgressBar progress("align");

        while (true) {
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

            const std::size_t chunk_size = chunk.size();

            // 刷新所有 writer
            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }

            std::vector<seq_io::SeqRecord>().swap(chunk);
            progress.tick(chunk_size);
        }

        // 完成并确保所有数据写入磁盘
        progress.done();
        spdlog::info("Alignment completed");

        for (auto& w : outs) {
            if (w) w->flush();
        }
        for (auto& w : outs_with_insertion) {
            if (w) w->flush();
        }
    }

    void RefAligner::alignSeq2Profile(const FilePath& qry_fasta_path, std::size_t batch_size)
    {
        profile_alignment_mode = true;

        // 参数检查和初始化
        if (ref_sequences.empty() || ref_sketch.empty()) {
            throw std::runtime_error("RefAligner::alignQueryToRef: reference sequence is empty");
        }

        constexpr std::size_t default_batch_size = 1000;
        if (batch_size == 0) {
            batch_size = default_batch_size;
        }

        // 设置线程数
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
        const int nthreads = std::max(1, omp_get_max_threads());

        const FilePath result_dir = work_dir / RESULTS_DIR;
        file_io::ensureDirectoryExists(result_dir, "result directory");

        spdlog::info("Starting alignment: {} threads, batch size {}", nthreads, batch_size);

        // 为每个线程创建独立的输出文件
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

        // 流式读取 + 批处理并行
        seq_io::KseqReader reader(qry_fasta_path);
        std::vector<seq_io::SeqRecord> chunk;
        std::vector<cigar::Cigar_t> cigar_chunk;
        std::vector<std::string> ref_id_chunk;
        chunk.reserve(batch_size);
        ProgressBar progress("align", 10);
        progress.tick(0);



        // 预热阶段：先串行处理固定数量的序列，并且“每比对一条就更新一次 profile”。
        // 目的：让后续 batch 并行阶段在更有信息量的 profile 上工作，降低冷启动阶段的偏差。
        // constexpr std::size_t profile_warmup_count = 5000;
        // std::vector<seq_io::SeqRecord> warmup_chunk(1);
        // std::vector<cigar::Cigar_t> warmup_cigar(1);
        // std::vector<std::string> warmup_ref_id(1);
        // spdlog::info("Warming up alignment with {} sequences, it may be slow", profile_warmup_count);
        //
        // std::size_t warmup_processed = 0;
        // seq_io::SeqRecord warmup_rec;
        // auto& warmup_out = *outs[0];
        // auto& warmup_out_insertion = *outs_with_insertion[0];
        //
        // while (warmup_processed < profile_warmup_count && reader.next(warmup_rec)) {
        //     warmup_chunk[0] = std::move(warmup_rec);
        //
        //     // 串行比对一条，得到该条最终使用的 CIGAR/参考索引。
        //     alignOneQueryToProfile(
        //         warmup_chunk[0],
        //         warmup_out,
        //         warmup_out_insertion,
        //         warmup_cigar[0],
        //         warmup_ref_id[0],
        //         threads);
        //
        //     // 每条序列比对完成后立即更新一次 profile，严格满足“比对一次、更新一次”。
        //     updateProfilesFromChunk(warmup_chunk, warmup_cigar, warmup_ref_id);
        //
        //     warmup_cigar[0].clear();
        //     warmup_ref_id[0].clear();
        //     ++warmup_processed;
        //     progress.tick();
        // 	if (warmup_processed % 10 == 0) {
        // 		spdlog::info("Warmup processed: {}/{} sequences", warmup_processed, profile_warmup_count);
        // 	}
        // }
        //
        // // 预热阶段统一刷新一次，避免仅 0 号 writer 长时间缓存。
        // warmup_out.flush();
        // warmup_out_insertion.flush();
        //
        // spdlog::info("alignSeq2Profile warmup processed {} sequences", warmup_processed);

        while (true) {
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

            // 为本批次预分配结果槽位，保证并行区内“每个 i 独占写入”无竞态
            cigar_chunk.clear();
            cigar_chunk.resize(chunk.size());
            ref_id_chunk.assign(chunk.size(), std::string());

#pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk, cigar_chunk, ref_id_chunk)
            {
                const int tid = omp_get_thread_num();
                auto& out = *outs[static_cast<std::size_t>(tid)];
                auto& out_insertion = *outs_with_insertion[static_cast<std::size_t>(tid)];

#pragma omp for schedule(dynamic, 1)

                for (std::int64_t i = 0; i < static_cast<std::int64_t>(chunk.size()); ++i) {
                    alignOneQueryToProfile(
                        chunk[static_cast<std::size_t>(i)],
                        out,
                        out_insertion,
                        cigar_chunk[static_cast<std::size_t>(i)],
                        ref_id_chunk[static_cast<std::size_t>(i)], 0);
                }
            }

            // 根据本批次对齐结果增量更新 profile，供后续批次选择参考和 profile 比对使用。
            updateProfilesFromChunk(chunk, cigar_chunk, ref_id_chunk);

            const std::size_t chunk_size = chunk.size();

            // 刷新所有 writer
            for (auto& w : outs) {
                w->flush();
            }
            for (auto& w : outs_with_insertion) {
                w->flush();
            }

            std::vector<seq_io::SeqRecord>().swap(chunk);
            progress.tick(chunk_size);
        }

        // 完成并确保所有数据写入磁盘
        progress.done();
        spdlog::info("Alignment completed");

        for (auto& w : outs) {
            if (w) w->flush();
        }
        for (auto& w : outs_with_insertion) {
            if (w) w->flush();
        }
    }


} // namespace align

