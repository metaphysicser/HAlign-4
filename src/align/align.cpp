#include "align.h"
#include "ksw2.h"

extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}

// 序列比对算法封装：KSW2 / WFA2
// - 统一返回 cigar::Cigar_t（压缩格式）
// - 支持全局比对、延伸比对、锚点分段比对

namespace align
{
    // KSW2 全局比对（end-to-end）- 编码序列并调用 KSW2
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query)
    {
        align::AlignConfig cfg;
        return globalAlignKSW2(ref, query, cfg);
    }

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query,
                                   align::AlignConfig cfg)
    {
        // 边界：任意一条序列为空时，直接返回纯 I / 纯 D 的 CIGAR
        if (ref.size() == 0 || query.size() == 0) {
            cigar::Cigar_t cigar;
            if (ref.size() == 0 && query.size() > 0) {
                cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(query.size())));
            } else if (query.size() == 0 && ref.size() > 0) {
                cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(ref.size())));
            }
            return cigar;
        }

        // 编码序列：DNA5 (A/C/G/T/N -> 0..4)
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());

        for (size_t i = 0; i < ref.size(); ++i)
            ref_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        cfg.band_width = align::auto_band(ref.size(), query.size());

        // 调用 KSW2
        ksw_extz_t ez{};
        ksw_extz2_sse(0,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        // 拷贝并释放 CIGAR
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);
        return cigar;
    }

    // KSW2 延伸比对 - 使用 zdrop 和 EXTZ_ONLY 标志
    cigar::Cigar_t extendAlignKSW2(const std::string& ref,
        const std::string& query,
        int zdrop)
    {
        // 编码序列
        std::vector<uint8_t> ref_enc(ref.size());
        std::vector<uint8_t> qry_enc(query.size());
        for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = align::ScoreChar2Idx[(uint8_t)ref[i]];
        for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = align::ScoreChar2Idx[(uint8_t)query[i]];

        // 配置参数：EXTZ_ONLY + RIGHT + APPROX_DROP for extension
        align::AlignConfig cfg;
        cfg.mat = align::dna5_simd_mat;
        cfg.zdrop = zdrop;
        cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
        cfg.end_bonus = 50;
        cfg.alphabet_size = 5;
        cfg.gap_open = 6;
        cfg.gap_extend = 2;
        cfg.band_width = align::auto_band(ref.size(), query.size());

        // 调用 KSW2
        ksw_extz_t ez{};
        ksw_extz2_sse(nullptr,
            static_cast<int>(qry_enc.size()), qry_enc.data(),
            static_cast<int>(ref_enc.size()), ref_enc.data(),
            cfg.alphabet_size, cfg.mat,
            cfg.gap_open, cfg.gap_extend,
            cfg.band_width, cfg.zdrop, cfg.end_bonus,
            cfg.flag, &ez);

        // 拷贝并释放 CIGAR
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(ez.cigar);
        return cigar;
    }


    // WFA2 全局比对 - 使用 gap_affine 模式
    cigar::Cigar_t globalAlignWFA2(const std::string& ref,
        const std::string& query)
    {
        // 构建 WFA2 属性
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.mismatch = 3;
        attributes.affine_penalties.gap_opening = 4;
        attributes.affine_penalties.gap_extension = 1;
        attributes.memory_mode = wavefront_memory_high;

        // 创建并执行 aligner
        wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
        wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());

        // 提取 CIGAR
        uint32_t* cigar_buffer = nullptr;
        int cigar_length = 0;
        cigar_get_CIGAR(wf_aligner->cigar, false, &cigar_buffer, &cigar_length);

        cigar::Cigar_t cigar;
        cigar.reserve(static_cast<std::size_t>(cigar_length));
        for (int i = 0; i < cigar_length; ++i)
            cigar.push_back(cigar_buffer[i]);

        wavefront_aligner_delete(wf_aligner);
        return cigar;
    }

    // ------------------------------------------------------------------
    // extendAlignWFA2：WFA2 延伸比对
    // ------------------------------------------------------------------
    // WFA2 的延伸比对（ends-free extension）：
    // - 通常用于从种子位置向外延伸，快速获得局部比对结果
    // - 可配置 zdrop 阈值，控制延伸过程中的提前终止
    //
    // 本实现使用 wavefront_aligner_attr_default 的高内存模式：
    // - 适合长序列比对，但可能导致较高的内存占用
    // ------------------------------------------------------------------
    // cigar::Cigar_t extendAlignWFA2(const std::string& ref,
    //     const std::string& query, int zdrop)
    // {
    //     wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    //     attributes.distance_metric = gap_affine;
    //     attributes.affine_penalties.mismatch = 2;      // X > 0
    //     attributes.affine_penalties.gap_opening = 3;   // O >= 0
    //     attributes.affine_penalties.gap_extension = 1; // E > 0
    //     attributes.memory_mode = wavefront_memory_high;
    //     attributes.heuristic.strategy = wf_heuristic_zdrop;
    //     attributes.heuristic.zdrop = zdrop;
    //     attributes.heuristic.steps_between_cutoffs = 1;
    //     //// Create a WFAligner
    //     //
    //     wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    //
    //     wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
    //     /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
    //
    //     aligner.alignEnd2End(ref, query);*/
    //
    //     uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
    //     int cigar_length = 0; // Length of the CIGAR string.
    //     // Retrieve the CIGAR string from the wavefront aligner.
    //     cigar_get_CIGAR(wf_aligner->cigar, true, &cigar_buffer, &cigar_length);
    //
    //     /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    //     cigar::Cigar_t cigar;
    //
    //     for (int i = 0; i < cigar_length; ++i)
    //         cigar.push_back(cigar_buffer[i]);
    //
    //     wavefront_aligner_delete(wf_aligner);
    //
    //     return cigar;
    // }

    cigar::Cigar_t globalAlignPSW(const ProfileMatrix& ref, const std::string& query, align::AlignConfig cfg)
    {
        // 边界：任意一条序列为空时，直接返回纯 I / 纯 D 的 CIGAR
        if (ref.len == 0 || query.size() == 0) {
            cigar::Cigar_t cigar;
            if (ref.len == 0 && query.size() > 0) {
                cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(query.size())));
            } else if (query.size() == 0 && ref.len > 0) {
                cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(ref.len)));
            }
            return cigar;
        }

        // 编码序列：DNA5 (A/C/G/T/N -> 0..4)

        std::vector<uint8_t> qry_enc(query.size());


        for (size_t i = 0; i < query.size(); ++i)
            qry_enc[i] = align::ScoreChar2Idx[static_cast<uint8_t>(query[i])];

        cfg.band_width = align::auto_band(ref.len, query.size());

        int m_cigar, n_cigar;
        uint32_t* cigar1 = nullptr;
        psw_prof_t ref_prof;
        ref_prof.len = ref.len;
        ref_prof.dim = ref.dim;
        ref_prof.depth = ref.depth;
        ref_prof.prof = ref.prof.data();


        float score = psw_gg3_sse_ps(0, query.size(), qry_enc.data(), ref.len, &ref_prof, (int8_t)ref.dim, cfg.mat,
                                    cfg.gap_open, cfg.gap_extend, cfg.band_width,
                                    &m_cigar, &n_cigar, &cigar1);

        // 拷贝并释放 CIGAR
        cigar::Cigar_t cigar;
        cigar.reserve(n_cigar);
        for (int i = 0; i < n_cigar; ++i)
            cigar.push_back(cigar1[i]);

        free(cigar1);
        return cigar;
    }

    // 基于锚点的分段全局比对（minimap2 风格）
    // - 用锚点拆分为多个片段，逐段全局比对，最后合并
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors)
    {
        align::AlignConfig cfg;
        align::AlignConfig first_cfg;
        first_cfg.flag = KSW_EZ_GENERIC_SC;

        const std::size_t ref_len = ref.size();
        const std::size_t qry_len = query.size();

        // 链化锚点：获取最佳链
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return globalAlignKSW2(ref, query);
        }

        // 按 query 坐标排序
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        cigar::Cigar_t result;
        result.reserve(chain_anchors.size() * 2 + 2);

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto append_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                  std::size_t qry_start, std::size_t qry_end, align::AlignConfig seg_cfg) {
            // 边界裁剪
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            const std::string seg_ref = ref.substr(ref_start, ref_end - ref_start);
            const std::string seg_qry = query.substr(qry_start, qry_end - qry_start);

            cigar::Cigar_t seg_cigar = globalAlignKSW2(seg_ref, seg_qry, seg_cfg);

            // 用 CIGAR 反推消耗长度
            const std::size_t seg_ref_len = seg_ref.size();
            const std::size_t seg_qry_len = seg_qry.size();
            const std::size_t c_ref = cigar::getRefLength(seg_cigar);
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar);

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignMM2(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing robust fallback for this segment",
                             seg_ref_len, seg_qry_len, c_ref, c_qry);
#endif
                // 兜底策略：Query 全 I、Ref 全 D
                cigar::Cigar_t forced_cigar;
                if (seg_qry_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry_len)));
                }
                if (seg_ref_len > 0) {
                    forced_cigar.push_back(cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref_len)));
                }
                cigar::appendCigar(result, forced_cigar);

                ref_pos = ref_end;
                qry_pos = qry_end;
                return;
            }

            cigar::appendCigar(result, seg_cigar);

            ref_pos = ref_start + c_ref;
            qry_pos = qry_start + c_qry;
        };

        // 左端：起点到第一个锚点
        {
            const auto& first = chain_anchors.front();
            append_segment(ref_pos, first.pos_ref, qry_pos, first.pos_qry, first_cfg);
        }

        // 逐锚点：处理 span 和 gap
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& a = chain_anchors[i];

            const std::size_t a_ref_start = static_cast<std::size_t>(a.pos_ref);
            const std::size_t a_qry_start = static_cast<std::size_t>(a.pos_qry);
            const std::size_t a_ref_end = a_ref_start + static_cast<std::size_t>(a.span);
            const std::size_t a_qry_end = a_qry_start + static_cast<std::size_t>(a.span);

            append_segment(ref_pos, a_ref_end, qry_pos, a_qry_end, cfg);

            if (i + 1 < chain_anchors.size()) {
                const auto& b = chain_anchors[i + 1];
                append_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry, cfg);
            }
        }

        // 右端：最后一个锚点到末尾
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg);

        // 最终一致性检查
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignMM2: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignKSW2(ref, query);
        }

        return result;
    }

} // namespace align
