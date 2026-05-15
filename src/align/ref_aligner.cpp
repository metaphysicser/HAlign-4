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
#include <utility>
#include <vector>
#include <omp.h>
#include <unordered_map>
#include <cstdio>

namespace align {
    // 读取参考序列，计算索引，生成共识序列
    RefAligner::RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                           int minimizer_size, int minimizer_window,
                           int sketch_size, int sketch_kmer_size, bool noncanonical,
                           int threads, std::string msa_tool,
                           bool enable_wfa,
                           const FilePath& reference_msa_path,
                           std::array<int8_t, 25> score_matrix_in,
                           int gap_open,
                           int gap_extend,
                           int min_profile_references,
                           int max_profile_references,
                           double min_profile_reference_similarity,
                           bool auto_strand)
        : work_dir(work_dir),
          minimizer_size(minimizer_size),
          minimizer_window(minimizer_window),
          sketch_size(sketch_size),
          sketch_kmer_size(sketch_kmer_size),
          noncanonical(noncanonical),
          threads(threads),
          msa_tool(std::move(msa_tool)),
          enable_wfa(enable_wfa),
          score_matrix(score_matrix_in),
          gap_open(gap_open),
          gap_extend(gap_extend),
          min_profile_references(min_profile_references),
          max_profile_references(max_profile_references),
          min_profile_reference_similarity(min_profile_reference_similarity),
          auto_strand(auto_strand)
    {
        loadReference(ref_fasta_path, reference_msa_path);
    }

    void RefAligner::loadReference(const FilePath& ref_fasta_path, const FilePath& reference_msa_path)
    {
        // 加载参考序列，sketch 在读取完成后并行构建，避免串行 I/O 循环承担重计算。
        seq_io::KseqReader reader(ref_fasta_path);
        seq_io::SeqRecord rec;
        spdlog::info("Loading reference sequences from {}", ref_fasta_path.string());

        std::vector<seq_io::SeqRecord> ref_records;
        while (reader.next(rec)) {
            ref_records.push_back(std::move(rec));
        }
        const int build_threads = std::max(1, threads > 0 ? threads : omp_get_max_threads());
        std::vector<mash::Sketch> ref_sketches(ref_records.size());

#pragma omp parallel for default(none) schedule(dynamic) num_threads(build_threads) \
    shared(ref_records, ref_sketches) firstprivate(sketch_kmer_size, sketch_size, noncanonical, random_seed)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(ref_records.size()); ++i) {
            ref_sketches[static_cast<std::size_t>(i)] = mash::sketchFromSequence(
                ref_records[static_cast<std::size_t>(i)].seq,
                static_cast<std::size_t>(sketch_kmer_size),
                static_cast<std::size_t>(sketch_size),
                noncanonical,
                random_seed);
        }

        ref_sequences.reserve(ref_records.size());
        ref_sketch.reserve(ref_records.size());
        for (std::size_t i = 0; i < ref_records.size(); ++i) {
            auto [seq_it, inserted] = ref_sequences.emplace(ref_records[i].id, std::move(ref_records[i]));
            if (inserted) {
                ref_sketch.emplace(seq_it->first, std::move(ref_sketches[i]));
            }
        }
        if (ref_sketch.size() > 1) {
            const std::size_t removed_common_hashes = mash::removeCommonHashesFromSketches(ref_sketch);
            if (removed_common_hashes != 0) {
                spdlog::info("Removed {} hashes shared by all reference sketches", removed_common_hashes);
            }
        }
        if (ref_sketch.size() > kSketchBitsetRefThreshold) {
            ref_sketch_bitset_index.build(ref_sketch, false);
            spdlog::info("Built reference sketch bitset index: {} refs, {} hash dictionary entries",
                         ref_sketch_bitset_index.size(),
                         ref_sketch_bitset_index.dictionarySize());
        }
        spdlog::info("Loaded {} reference sequences", ref_sequences.size());

        // 设置共识序列生成的文件路径
        const bool has_prealigned_ref = !reference_msa_path.empty();
        const FilePath consensus_unaligned_file = has_prealigned_ref ? reference_msa_path : ref_fasta_path;
        const FilePath consensus_aligned_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath consensus_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
        const FilePath consensus_json_file = FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

        // 执行 MSA 并生成共识序列
        // 说明：
        // - 默认路径：对 -r/--reference 提供的原始 FASTA 重新做一次 MSA，再从对齐结果生成共识；
        // - 若用户同时提供 --reference-msa，则说明这个参考集合已经有现成的 MSA，
        //   这里直接复用该 MSA，避免重复对齐，从而节省时间并保持参考坐标不变。
        constexpr std::size_t consensus_batch_size = 4096;
        if (has_prealigned_ref) {
            file_io::copyFile(consensus_unaligned_file, consensus_aligned_file);
        } else {
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, this->msa_tool, threads);
        }

        seq_io::KseqReader reader2(consensus_aligned_file);
        seq_io::SeqRecord rec2;
        std::vector<seq_io::SeqRecord> aligned_ref_records;
        while (reader2.next(rec2)) {
            aligned_ref_records.push_back(std::move(rec2));
        }
        std::vector<ProfileMatrix> aligned_ref_profiles(aligned_ref_records.size());

#pragma omp parallel for default(none) schedule(dynamic) num_threads(build_threads) \
    shared(aligned_ref_records, aligned_ref_profiles)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(aligned_ref_records.size()); ++i) {
            aligned_ref_profiles[static_cast<std::size_t>(i)] =
                ProfileMatrix::fromAlignedSequence(aligned_ref_records[static_cast<std::size_t>(i)].seq);
        }

        ref_profile.reserve(aligned_ref_records.size());
        for (std::size_t i = 0; i < aligned_ref_records.size(); ++i) {
            ref_profile.emplace(aligned_ref_records[i].id, std::move(aligned_ref_profiles[i]));
        }
        rebuildNormalizedReferenceProfiles();
        consensus::ConsensusResult consensus_result = consensus::generateConsensusResult(
            consensus_aligned_file, consensus_file, consensus_json_file,
            0, threads, consensus_batch_size);

        consensus_gap_seq.id = "consensus";
        consensus_gap_seq.seq = std::move(consensus_result.gap_seq);
        consensus_profile = ProfileMatrix::fromConsensusCounts(consensus_result.counts);

        // 以下变量为 seq2seq 准备。
        consensus_seq.id = "consensus";
        consensus_seq.seq = std::move(consensus_result.seq);
        // 预计算共识序列的 sketch 和 minimizer，避免重复计算。
        consensus_sketch = mash::sketchFromSequence(
            consensus_seq.seq,
            static_cast<std::size_t>(sketch_kmer_size),
            static_cast<std::size_t>(sketch_size),
            noncanonical,
            random_seed);

        consensus_minimizer = minimizer::extractMinimizer(
            consensus_seq.seq, minimizer_size, minimizer_window, noncanonical);
        consensus_gap_minimizer = minimizer::extractMinimizer(
            consensus_gap_seq.seq, minimizer_size, minimizer_window, noncanonical);
    }

    // Options 配置委托构造
    RefAligner::RefAligner(const Options& opt, const FilePath& ref_fasta_path)
        : RefAligner(
            opt.workdir,
            ref_fasta_path,
            opt.minimizer_size,
            opt.minimizer_window,
            opt.sketch_size,
            opt.sketch_kmer_size,
            true,
            opt.threads,
            opt.msa_tool,
            opt.enable_wfa,
            FilePath(opt.reference_msa_path),
            opt.score_matrix,
            opt.gap_open,
            opt.gap_extend,
            opt.min_profile_references,
            opt.max_profile_references,
            opt.min_profile_reference_similarity,
            opt.auto_strand)
    {
    }

    void RefAligner::rebuildNormalizedReferenceProfiles()
    {
        std::vector<std::pair<std::string, const ProfileMatrix*>> profiles;
        profiles.reserve(ref_profile.size());
        for (const auto& [id, profile] : ref_profile) {
            profiles.emplace_back(id, &profile);
        }

        const int build_threads = std::max(1, threads > 0 ? threads : omp_get_max_threads());
        std::vector<ProfileMatrix> normalized_profiles(profiles.size());

#pragma omp parallel for default(none) schedule(dynamic) num_threads(build_threads) shared(profiles, normalized_profiles)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(profiles.size()); ++i) {
            normalized_profiles[static_cast<std::size_t>(i)] =
                normalizeProfileEqualWeight(*profiles[static_cast<std::size_t>(i)].second);
        }

        normalized_ref_profile.clear();
        normalized_ref_profile.reserve(profiles.size());
        for (std::size_t i = 0; i < profiles.size(); ++i) {
            normalized_ref_profile.emplace(profiles[i].first, std::move(normalized_profiles[i]));
        }
    }

    void RefAligner::refreshNormalizedReferenceProfiles(const std::vector<std::string>& ref_ids)
{
    if (ref_ids.empty()) {
        return;
    }

    std::vector<std::string> unique_ids = ref_ids;
    std::sort(unique_ids.begin(), unique_ids.end());
    unique_ids.erase(std::unique(unique_ids.begin(), unique_ids.end()), unique_ids.end());

    // 避免后面 emplace 时 rehash。
    normalized_ref_profile.reserve(ref_profile.size());

    // 先确保 normalized_ref_profile 中存在对应 key。
    // 注意：不要在并行区里访问 / 修改 unordered_map。
    for (const std::string& ref_id : unique_ids) {
        if (normalized_ref_profile.find(ref_id) == normalized_ref_profile.end()) {
            normalized_ref_profile.emplace(ref_id, ProfileMatrix());
        }
    }

    struct NormalizeTask {
        const ProfileMatrix* src = nullptr;
        ProfileMatrix* dst = nullptr;
    };

    std::vector<NormalizeTask> tasks;
    tasks.reserve(unique_ids.size());

    for (const std::string& ref_id : unique_ids) {
        auto src_it = ref_profile.find(ref_id);
        if (src_it == ref_profile.end()) {
            throw std::runtime_error("Reference profile '" + ref_id + "' not found");
        }

        auto dst_it = normalized_ref_profile.find(ref_id);
        if (dst_it == normalized_ref_profile.end()) {
            throw std::runtime_error("Normalized reference profile '" + ref_id + "' not found");
        }

        const ProfileMatrix& src = src_it->second;

        if (src.len < 0 || src.dim <= 0) {
            throw std::runtime_error("refreshNormalizedReferenceProfiles: invalid profile shape for '" + ref_id + "'");
        }

        const std::size_t n_cells =
            static_cast<std::size_t>(src.len) * static_cast<std::size_t>(src.dim);

        if (src.prof.size() != n_cells) {
            throw std::runtime_error("refreshNormalizedReferenceProfiles: profile size mismatch for '" + ref_id + "'");
        }

        tasks.push_back(NormalizeTask{&src_it->second, &dst_it->second});
    }

    auto is_power_of_two = [](std::uint64_t x) -> bool {
        return x != 0 && (x & (x - 1)) == 0;
    };

    auto log2_u64 = [](std::uint64_t x) -> unsigned {
        unsigned shift = 0;
        while (x > 1) {
            x >>= 1;
            ++shift;
        }
        return shift;
    };

    auto normalize_into = [is_power_of_two, log2_u64](
        const ProfileMatrix& src,
        ProfileMatrix& dst)
    {
        const std::size_t n_cells =
            static_cast<std::size_t>(src.len) * static_cast<std::size_t>(src.dim);

        dst.len = src.len;
        dst.dim = src.dim;
        dst.depth = static_cast<int>(kProfileEqualWeightScale);

        // 关键：resize 会复用已有 capacity，避免每次 refresh 都重新分配大块内存。
        dst.prof.resize(n_cells);

        const std::uint32_t* src_data = src.prof.data();
        std::uint32_t* dst_data = dst.prof.data();

        const std::uint64_t depth =
            src.depth > 0 ? static_cast<std::uint64_t>(src.depth) : 1ULL;

        if (depth == static_cast<std::uint64_t>(kProfileEqualWeightScale)) {
            std::copy(src_data, src_data + n_cells, dst_data);
            return;
        }

        if (is_power_of_two(depth)) {
            const unsigned depth_shift = log2_u64(depth);

            for (std::size_t i = 0; i < n_cells; ++i) {
                std::uint64_t numerator;

                if constexpr (kProfileEqualWeightScale == 1024) {
                    numerator = static_cast<std::uint64_t>(src_data[i]) << 10;
                } else {
                    numerator = static_cast<std::uint64_t>(src_data[i]) *
                                static_cast<std::uint64_t>(kProfileEqualWeightScale);
                }

                numerator += depth >> 1;

                dst_data[i] = static_cast<std::uint32_t>(numerator >> depth_shift);
            }

            return;
        }

        for (std::size_t i = 0; i < n_cells; ++i) {
            std::uint64_t numerator =
                static_cast<std::uint64_t>(src_data[i]) *
                static_cast<std::uint64_t>(kProfileEqualWeightScale);

            numerator += depth / 2;

            dst_data[i] = static_cast<std::uint32_t>(numerator / depth);
        }
    };

    const int build_threads = std::max(1, threads > 0 ? threads : omp_get_max_threads());

#pragma omp parallel for default(none) schedule(static) num_threads(build_threads) \
    shared(tasks, normalize_into)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(tasks.size()); ++i) {
        NormalizeTask& task = tasks[static_cast<std::size_t>(i)];
        normalize_into(*task.src, *task.dst);
    }
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
            ref_mz_tmp = minimizer::extractMinimizer(ref, minimizer_size, minimizer_window, noncanonical);
            ref_mz_ptr = &ref_mz_tmp;
        }
        if (qry_mz_ptr == nullptr || qry_mz_ptr->empty()) {
            qry_mz_tmp = minimizer::extractMinimizer(query, minimizer_size, minimizer_window, noncanonical);
            qry_mz_ptr = &qry_mz_tmp;
        }

        const anchor::Anchors anchors = minimizer::collect_anchors(*ref_mz_ptr, *qry_mz_ptr);
        cigar::Cigar_t result = enable_wfa
            ? globalAlignWFA2(ref, query)
            : globalAlignSeq2Seq(ref, query, anchors, makeAlignConfig());

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
            ref_mz_tmp = minimizer::extractMinimizer(ref_string, minimizer_size, minimizer_window, noncanonical);
            ref_mz_ptr = &ref_mz_tmp;
        }
        if (qry_mz_ptr == nullptr || qry_mz_ptr->empty()) {
            qry_mz_tmp = minimizer::extractMinimizer(query, minimizer_size, minimizer_window, noncanonical);
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

    std::vector<mash::SketchMatch> RefAligner::selectProfileMatches(const mash::Sketch& query_sketch) const
    {
        std::vector<mash::SketchMatch> matches;
        if (ref_sketch_bitset_index.usable()) {
            matches = ref_sketch_bitset_index.findAll(query_sketch);
        } else {
            matches.reserve(ref_sketch.size());
            for (const auto& [ref_id, sketch] : ref_sketch) {
                matches.push_back(mash::SketchMatch{
                    ref_id,
                    mash::jaccard(query_sketch, sketch),
                    true
                });
            }
        }

        const auto clamp01 = [](double value) noexcept {
            if (value < 0.0) return 0.0;
            if (value > 1.0) return 1.0;
            return value;
        };

        const std::size_t min_count = static_cast<std::size_t>(std::max(1, min_profile_references));
        const std::size_t max_count = std::max(
            min_count,
            static_cast<std::size_t>(std::max(min_profile_references, max_profile_references)));
        const double min_sequence_similarity = clamp01(min_profile_reference_similarity);
        const std::size_t mash_kmer_len = static_cast<std::size_t>(std::max(1, sketch_kmer_size));

        // for (mash::SketchMatch& match : matches) {
        //     match.similarity = mash::aniFromJaccard(clamp01(match.similarity), mash_kmer_len);
        // }

        std::sort(matches.begin(), matches.end(),
                  [](const mash::SketchMatch& a, const mash::SketchMatch& b) {
                      if (a.similarity != b.similarity) return a.similarity > b.similarity;
                      return a.id < b.id;
                  });

        std::vector<mash::SketchMatch> selected;
        selected.reserve(std::min(matches.size(), max_count));

        for (mash::SketchMatch& match : matches) {
            if (selected.size() < min_count) {
                selected.push_back(std::move(match));
                continue;
            }

            if (selected.size() >= max_count ||
                match.similarity < min_sequence_similarity) {
                break;
            }

            selected.push_back(std::move(match));
        }

        return selected;
    }

    // 单条 query 比对
    void RefAligner::alignOneQueryToRef(const seq_io::SeqRecord& q,
                                       seq_io::SeqWriter& out,
                                       seq_io::SeqWriter& out_insertion) const
    {
        const seq_io::SeqRecord* query_record = &q;
        seq_io::SeqRecord rc_query_storage;
        if (auto_strand) {
            const mash::Sketch fwd_sketch = mash::sketchFromSequence(
                q.seq,
                static_cast<std::size_t>(sketch_kmer_size),
                static_cast<std::size_t>(sketch_size),
                noncanonical,
                random_seed);
            rc_query_storage = reverseComplementRecord(q);
            const mash::Sketch rc_sketch = mash::sketchFromSequence(
                rc_query_storage.seq,
                static_cast<std::size_t>(sketch_kmer_size),
                static_cast<std::size_t>(sketch_size),
                noncanonical,
                random_seed);
            if (mash::jaccard(rc_sketch, consensus_sketch) >
                mash::jaccard(fwd_sketch, consensus_sketch)) {
                query_record = &rc_query_storage;
            }
        }

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            query_record->seq, minimizer_size, minimizer_window, noncanonical);

        const double consensus_similarity = 1;
        cigar::Cigar_t consensus_cigar = Seq2SeqWithAnchor(
            consensus_seq.seq, query_record->seq, consensus_similarity,
            &consensus_minimizer, &query_minimizer);

        if (cigar::hasInsertion(consensus_cigar)) {
            writeSamRecord(*query_record, consensus_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(*query_record, consensus_cigar, consensus_seq.id, out);
        }
    }

    void RefAligner::alignOneQueryToProfile(const seq_io::SeqRecord& q,
                       seq_io::SeqWriter& out,
                       seq_io::SeqWriter& out_insertion,
                       cigar::Cigar_t& out_cigar,
                       std::vector<std::string>& out_ref_ids,
                       seq_io::SeqRecord& out_profile_query,
                       int thread) const
    {
        // 初始化输出占位，便于调用方在调试时识别“未写入”状态
        out_cigar.clear();
        out_ref_ids.clear();
        out_profile_query = q;

        auto sketch_query = [this](const std::string& seq) {
            return mash::sketchFromSequence(
                seq,
                static_cast<std::size_t>(sketch_kmer_size),
                static_cast<std::size_t>(sketch_size),
                noncanonical,
                random_seed);
        };

        mash::Sketch qsk = sketch_query(q.seq);
        std::vector<mash::SketchMatch> selected_matches;
        double best_sequence_similarity = -1.0;

        if (ref_sequences.size() > 1) {
            selected_matches = selectProfileMatches(qsk);
            best_sequence_similarity = selected_matches.empty() ? 0.0 : selected_matches.front().similarity;

            if (auto_strand) {
                seq_io::SeqRecord rc_query = reverseComplementRecord(q);
                mash::Sketch rc_sketch = sketch_query(rc_query.seq);
                std::vector<mash::SketchMatch> rc_matches = selectProfileMatches(rc_sketch);
                const double rc_best = rc_matches.empty() ? 0.0 : rc_matches.front().similarity;
                if (rc_best > best_sequence_similarity) {
                    out_profile_query = std::move(rc_query);
                    qsk = std::move(rc_sketch);
                    selected_matches = std::move(rc_matches);
                    best_sequence_similarity = rc_best;
                }
            }
        } else {
            best_sequence_similarity = mash::aniFromJaccard(
                mash::jaccard(qsk, consensus_sketch),
                static_cast<std::size_t>(std::max(1, sketch_kmer_size)));
            if (auto_strand) {
                seq_io::SeqRecord rc_query = reverseComplementRecord(q);
                mash::Sketch rc_sketch = sketch_query(rc_query.seq);
                const double rc_sequence_similarity = mash::aniFromJaccard(
                    mash::jaccard(rc_sketch, consensus_sketch),
                    static_cast<std::size_t>(std::max(1, sketch_kmer_size)));
                if (rc_sequence_similarity > best_sequence_similarity) {
                    out_profile_query = std::move(rc_query);
                    qsk = std::move(rc_sketch);
                    best_sequence_similarity = rc_sequence_similarity;
                }
            }
        }

        const SeedHits query_minimizer = minimizer::extractMinimizer(
            out_profile_query.seq, minimizer_size, minimizer_window, noncanonical);

        const ProfileMatrix* alignment_profile = &consensus_profile;
        std::string alignment_ref_string = consensus_gap_seq.seq;
        const SeedHits* alignment_ref_minimizer = &consensus_gap_minimizer;
        SeedHits alignment_ref_minimizer_storage;
        ProfileMatrix combined_profile;

        if (ref_sequences.size() > 1) {
            if (selected_matches.empty()) {
                throw std::runtime_error("alignOneQueryToProfile: no reference profile candidate selected");
            }

            std::vector<const ProfileMatrix*> selected_profiles;
            selected_profiles.reserve(selected_matches.size());
            out_ref_ids.reserve(selected_matches.size());

            for (const mash::SketchMatch& match : selected_matches) {
                auto profile_it = normalized_ref_profile.find(match.id);
                if (profile_it == normalized_ref_profile.end()) {
                    throw std::runtime_error("Normalized reference profile '" + match.id + "' not found");
                }
                selected_profiles.push_back(&profile_it->second);
                out_ref_ids.push_back(match.id);
            }

            combined_profile = combineProfilesEqualWeight(selected_profiles);
            alignment_profile = &combined_profile;
            alignment_ref_string = profileToGappedSequence(combined_profile);
            alignment_ref_minimizer_storage = minimizer::extractMinimizer(
                alignment_ref_string, minimizer_size, minimizer_window, noncanonical);
            alignment_ref_minimizer = &alignment_ref_minimizer_storage;
        }

        // 执行全局比对
        cigar::Cigar_t initial_cigar = Seq2ProfileWithAnchor(
            *alignment_profile, alignment_ref_string, out_profile_query.seq, best_sequence_similarity, thread,
            alignment_ref_minimizer, &query_minimizer);

        out_cigar = initial_cigar;

        if (cigar::hasInsertion(initial_cigar)) {
            writeSamRecord(out_profile_query, initial_cigar, consensus_seq.id, out_insertion);
        } else {
            writeSamRecord(out_profile_query, initial_cigar, consensus_seq.id, out);
        }
    }

    bool RefAligner::applyCigarToProfile(
    const std::string& query_seq,
    const cigar::Cigar_t& cigar,
    ProfileMatrix& target_profile)
	{
	    // 关键约束：这里只做“计数累加”，不改变 profile 的列数/维度，避免影响现有比对逻辑。
	    if (target_profile.len < 0 || target_profile.dim <= 0) {
	        return false;
	    }

	    const std::size_t profile_len = static_cast<std::size_t>(target_profile.len);
	    const std::size_t profile_dim = static_cast<std::size_t>(target_profile.dim);
	    const std::size_t query_len = query_seq.size();

	    if (profile_dim == 0 || target_profile.prof.size() != profile_len * profile_dim) {
	        return false;
	    }

	    std::uint32_t* prof_data = target_profile.prof.data();
	    const unsigned char* query_data =
	        reinterpret_cast<const unsigned char*>(query_seq.data());

	    std::size_t ref_pos = 0;
	    std::size_t ref_offset = 0; // 等价于 ref_pos * profile_dim，但用递增避免循环内乘法
	    std::size_t qry_pos = 0;

	    for (const cigar::CigarUnit unit : cigar) {
	        char op = '\0';
	        std::uint32_t len = 0;
	        cigar::intToCigar(unit, op, len);

	        const std::size_t block_len = static_cast<std::size_t>(len);

	        switch (op) {
	            case 'M':
	            case '=':
	            case 'X': {
	                // M/=/X 同时消耗参考和 query。
	                // 原来是在每个碱基里检查 ref_pos/qry_pos 是否越界；
	                // 这里改成每个 CIGAR block 只检查一次。
	                if (ref_pos > profile_len ||
	                    qry_pos > query_len ||
	                    block_len > profile_len - ref_pos ||
	                    block_len > query_len - qry_pos) {
	                    return false;
	                }

	                for (std::size_t k = 0; k < block_len; ++k) {
	                    const std::size_t base_idx = static_cast<std::size_t>(
	                        align::ScoreChar2Idx[query_data[qry_pos]]);

	                    ++prof_data[ref_offset + base_idx];

	                    ref_offset += profile_dim;
	                    ++ref_pos;
	                    ++qry_pos;
	                }

	                break;
	            }

	            case 'D':
	            case 'N': {
	                // D/N 只消耗参考：该列没有 query 碱基贡献，仅推进参考坐标。
	                if (ref_pos > profile_len ||
	                    block_len > profile_len - ref_pos) {
	                    return false;
	                }

	                ref_pos += block_len;
	                ref_offset += block_len * profile_dim;

	                break;
	            }

	            case 'I':
	            case 'S': {
	                // I/S 只消耗 query：不对应参考列，不能写入 profile，直接推进 query 坐标。
	                if (qry_pos > query_len ||
	                    block_len > query_len - qry_pos) {
	                    return false;
	                }

	                qry_pos += block_len;

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
	    const std::vector<std::vector<std::string>>& ref_ids_chunk)
	{

	    std::size_t empty_cigar_count = 0;
	    std::size_t invalid_cigar_count = 0;
	    std::size_t consensus_update_count = 0;
	    std::size_t consensus_skip_count = 0;
	    std::size_t ref_not_found_count = 0;
	    std::size_t ref_update_count = 0;
	    std::size_t unique_dirty_ref_count = 0;
	    std::size_t total_ref_edges = 0;

	    {
	        if (chunk.size() != cigar_chunk.size() || chunk.size() != ref_ids_chunk.size()) {

	            throw std::runtime_error("updateProfilesFromChunk: chunk/cigar/ref_ids size mismatch");
	        }
	    }

	    std::size_t ref_profile_len = 0;
	    std::size_t ref_profile_dim = 0;
	    std::size_t consensus_profile_len = 0;
	    std::size_t consensus_profile_dim = 0;

	    {


	        if (!ref_profile.empty()) {
	            const ProfileMatrix& first_ref_profile = ref_profile.begin()->second;

	            if (first_ref_profile.len < 0 || first_ref_profile.dim <= 0) {
	                throw std::runtime_error("updateProfilesFromChunk: invalid reference profile shape");
	            }

	            ref_profile_len = static_cast<std::size_t>(first_ref_profile.len);
	            ref_profile_dim = static_cast<std::size_t>(first_ref_profile.dim);
	        }

	        if (consensus_profile.len >= 0 && consensus_profile.dim > 0) {
	            consensus_profile_len = static_cast<std::size_t>(consensus_profile.len);
	            consensus_profile_dim = static_cast<std::size_t>(consensus_profile.dim);
	        }

	    }

	    struct QueryUpdateOffsets {
	        bool valid = false;
	        bool consensus_target = false;
	        std::vector<std::size_t> offsets;
	    };

	    auto build_update_offsets = [](
	        const std::string& query_seq,
	        const cigar::Cigar_t& cigar,
	        std::size_t profile_len,
	        std::size_t profile_dim,
	        std::vector<std::size_t>& offsets) -> bool
	    {
	        if (profile_dim == 0) {
	            return false;
	        }

	        offsets.clear();
	        offsets.reserve(std::min(query_seq.size(), profile_len));

	        const std::size_t query_len = query_seq.size();
	        const unsigned char* query_data =
	            reinterpret_cast<const unsigned char*>(query_seq.data());

	        std::size_t ref_pos = 0;
	        std::size_t ref_offset = 0;
	        std::size_t qry_pos = 0;

	        for (const cigar::CigarUnit unit : cigar) {
	            char op = '\0';
	            std::uint32_t len = 0;
	            cigar::intToCigar(unit, op, len);

	            const std::size_t block_len = static_cast<std::size_t>(len);

	            switch (op) {
	                case 'M':
	                case '=':
	                case 'X': {
	                    if (ref_pos > profile_len ||
	                        qry_pos > query_len ||
	                        block_len > profile_len - ref_pos ||
	                        block_len > query_len - qry_pos) {
	                        return false;
	                    }

	                    for (std::size_t k = 0; k < block_len; ++k) {
	                        const std::size_t base_idx = static_cast<std::size_t>(
	                            align::ScoreChar2Idx[query_data[qry_pos]]);

	                        offsets.push_back(ref_offset + base_idx);

	                        ref_offset += profile_dim;
	                        ++ref_pos;
	                        ++qry_pos;
	                    }

	                    break;
	                }

	                case 'D':
	                case 'N': {
	                    if (ref_pos > profile_len ||
	                        block_len > profile_len - ref_pos) {
	                        return false;
	                    }

	                    ref_pos += block_len;
	                    ref_offset += block_len * profile_dim;

	                    break;
	                }

	                case 'I':
	                case 'S': {
	                    if (qry_pos > query_len ||
	                        block_len > query_len - qry_pos) {
	                        return false;
	                    }

	                    qry_pos += block_len;

	                    break;
	                }

	                case 'H':
	                case 'P': {
	                    break;
	                }

	                default:
	                    return false;
	            }
	        }

	        return ref_pos == profile_len;
	    };

	    auto apply_offsets_to_profile = [](
	        const std::vector<std::size_t>& offsets,
	        ProfileMatrix& profile)
	    {
	        std::uint32_t* prof_data = profile.prof.data();

	        for (const std::size_t offset : offsets) {
	            ++prof_data[offset];
	        }
	    };

	    std::vector<QueryUpdateOffsets> query_offsets(chunk.size());

	    {

	        const int update_threads = std::max(1, threads > 0 ? threads : omp_get_max_threads());

	#pragma omp parallel for default(none) schedule(static) num_threads(update_threads) \
	    shared(chunk, cigar_chunk, ref_ids_chunk, query_offsets, build_update_offsets) \
	    firstprivate(ref_profile_len, ref_profile_dim, consensus_profile_len, consensus_profile_dim)
	        for (std::int64_t si = 0; si < static_cast<std::int64_t>(chunk.size()); ++si) {
	            const std::size_t i = static_cast<std::size_t>(si);

	            if (cigar_chunk[i].empty()) {
	                continue;
	            }

	            QueryUpdateOffsets& qoff = query_offsets[i];

	            if (ref_ids_chunk[i].empty()) {
	                qoff.consensus_target = true;
	                qoff.valid = build_update_offsets(
	                    chunk[i].seq,
	                    cigar_chunk[i],
	                    consensus_profile_len,
	                    consensus_profile_dim,
	                    qoff.offsets);
	            } else {
	                qoff.consensus_target = false;
	                qoff.valid = build_update_offsets(
	                    chunk[i].seq,
	                    cigar_chunk[i],
	                    ref_profile_len,
	                    ref_profile_dim,
	                    qoff.offsets);
	            }
	        }

	        for (std::size_t i = 0; i < chunk.size(); ++i) {
	            if (cigar_chunk[i].empty()) {
	                ++empty_cigar_count;
	            } else if (!query_offsets[i].valid) {
	                ++invalid_cigar_count;
	            }
	        }

	    }

	    struct RefUpdateGroup {
	        std::string ref_id;
	        ProfileMatrix* profile = nullptr;
	        std::vector<std::size_t> query_indices;
	    };

	    std::vector<RefUpdateGroup> ref_groups;
	    std::unordered_map<std::string, std::size_t> ref_group_index;
	    std::vector<std::size_t> consensus_query_indices;

	    {

	        ref_groups.reserve(ref_profile.size());
	        ref_group_index.reserve(ref_profile.size());
	        consensus_query_indices.reserve(chunk.size());

	        for (std::size_t i = 0; i < chunk.size(); ++i) {
	            if (cigar_chunk[i].empty() || !query_offsets[i].valid) {
	                continue;
	            }

	            if (ref_ids_chunk[i].empty()) {
	                consensus_query_indices.push_back(i);
	                continue;
	            }

	            for (const std::string& ref_id : ref_ids_chunk[i]) {
	                ++total_ref_edges;

	                auto group_it = ref_group_index.find(ref_id);
	                if (group_it == ref_group_index.end()) {
	                    auto profile_it = ref_profile.find(ref_id);
	                    if (profile_it == ref_profile.end()) {
	                        ++ref_not_found_count;

	#ifdef _DEBUG
	                        spdlog::debug("updateProfilesFromChunk: skip invalid ref_id={} at i={}",
	                                      ref_id, i);
	#endif
	                        continue;
	                    }

	                    const std::size_t group_index = ref_groups.size();

	                    ref_group_index.emplace(ref_id, group_index);

	                    RefUpdateGroup group;
	                    group.ref_id = ref_id;
	                    group.profile = &profile_it->second;
	                    group.query_indices.reserve(64);
	                    group.query_indices.push_back(i);

	                    ref_groups.push_back(std::move(group));
	                } else {
	                    ref_groups[group_it->second].query_indices.push_back(i);
	                }
	            }
	        }

	        unique_dirty_ref_count = ref_groups.size();
	    }

	    {
	        for (const std::size_t query_index : consensus_query_indices) {
	            if (consensus_profile.len < 0 ||
	                consensus_profile.dim <= 0 ||
	                consensus_profile.prof.size() !=
	                    static_cast<std::size_t>(consensus_profile.len) *
	                    static_cast<std::size_t>(consensus_profile.dim)) {
	                ++consensus_skip_count;
	                continue;
	            }

	            apply_offsets_to_profile(query_offsets[query_index].offsets, consensus_profile);
	            ++consensus_profile.depth;
	            ++consensus_update_count;
	        }
	    }

		    {
    		const int update_threads = std::max(1, threads > 0 ? threads : omp_get_max_threads());

	#pragma omp parallel for default(none) schedule(dynamic, 1) num_threads(update_threads) \
	shared(ref_groups, query_offsets, apply_offsets_to_profile) \
	firstprivate(ref_profile_len, ref_profile_dim) \
	reduction(+:ref_update_count)
    		for (std::int64_t gi = 0; gi < static_cast<std::int64_t>(ref_groups.size()); ++gi) {
    			RefUpdateGroup& group = ref_groups[static_cast<std::size_t>(gi)];
    			ProfileMatrix& profile = *group.profile;

    			if (profile.len < 0 ||
					profile.dim <= 0 ||
					static_cast<std::size_t>(profile.len) != ref_profile_len ||
					static_cast<std::size_t>(profile.dim) != ref_profile_dim ||
					profile.prof.size() != ref_profile_len * ref_profile_dim) {
    				continue;
					}

    			for (const std::size_t query_index : group.query_indices) {
    				apply_offsets_to_profile(query_offsets[query_index].offsets, profile);
    			}

    			profile.depth += static_cast<int>(group.query_indices.size());
    			ref_update_count += group.query_indices.size();
    		}
		    }

	    std::vector<std::string> dirty_ref_ids;
	    dirty_ref_ids.reserve(ref_groups.size());
	    for (const RefUpdateGroup& group : ref_groups) {
	        dirty_ref_ids.push_back(group.ref_id);
	    }

	    {
	        refreshNormalizedReferenceProfiles(dirty_ref_ids);
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
        std::vector<seq_io::SeqRecord> profile_query_chunk;
        std::vector<cigar::Cigar_t> cigar_chunk;
        std::vector<std::vector<std::string>> ref_ids_chunk;
        chunk.reserve(batch_size);
        ProgressBar progress("align", 10);
        progress.tick(0);

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
            profile_query_chunk.clear();
            profile_query_chunk.resize(chunk.size());
            cigar_chunk.clear();
            cigar_chunk.resize(chunk.size());
            ref_ids_chunk.clear();
            ref_ids_chunk.resize(chunk.size());

#pragma omp parallel default(none) shared(outs, outs_with_insertion, chunk, profile_query_chunk, cigar_chunk, ref_ids_chunk)
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
                        ref_ids_chunk[static_cast<std::size_t>(i)],
                        profile_query_chunk[static_cast<std::size_t>(i)], 0);
                }
            }

            // 根据本批次对齐结果增量更新 profile，供后续批次选择参考和 profile 比对使用。
            updateProfilesFromChunk(profile_query_chunk, cigar_chunk, ref_ids_chunk);

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
