#include "align.h"
#include "ksw2.h"
#ifdef _OPENMP
#include <omp.h>
#endif

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
        //cfg.band_width = -1;

        // 改为调用 ksw_gg2_sse：该接口是标准全局比对（Needleman-Wunsch），
        // 直接返回完整路径 CIGAR，不再依赖 extz 的 zdrop/end_bonus/flag 行为。
        int m_cigar = 0;
        int n_cigar = 0;
        uint32_t* cigar_raw = nullptr;
        // ksw_gg2_sse(nullptr,
        //             static_cast<int>(qry_enc.size()), qry_enc.data(),
        //             static_cast<int>(ref_enc.size()), ref_enc.data(),
        //             static_cast<int8_t>(cfg.alphabet_size), cfg.mat,
        //             static_cast<int8_t>(cfg.gap_open), static_cast<int8_t>(cfg.gap_extend),
        //             cfg.band_width,
        //             &m_cigar, &n_cigar, &cigar_raw);
        ksw_extz_t ez{};
        ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
         cfg.flag, &ez);

        // 拷贝并释放 CIGAR：保持对外返回类型不变，避免调用方感知底层算法替换。
        // cigar::Cigar_t cigar;
        // if (n_cigar > 0 && cigar_raw != nullptr) {
        //     cigar.reserve(static_cast<std::size_t>(n_cigar));
        //     for (int i = 0; i < n_cigar; ++i)
        //         cigar.push_back(cigar_raw[i]);
        // }
        cigar::Cigar_t cigar;
        cigar.reserve(ez.n_cigar);
        for (int i = 0; i < ez.n_cigar; ++i)
            cigar.push_back(ez.cigar[i]);

        free(cigar_raw);
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
        //cfg.band_width = -1;

        int m_cigar = 0, n_cigar = 0;
        uint32_t* cigar1 = 0;
        psw_prof_t ref_prof;
        ref_prof.len = ref.len;
        ref_prof.dim = ref.dim;
        ref_prof.depth = ref.depth;
        ref_prof.prof = ref.prof.data();


        psw_gg3_sse_ps(0, query.size(), qry_enc.data(), ref.len, &ref_prof, (int8_t)ref.dim, cfg.mat,
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
    cigar::Cigar_t globalAlignSeq2Seq(const std::string& ref,
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
                spdlog::warn("globalAlignSeq2Seq(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing robust fallback for this segment",
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
            spdlog::error("globalAlignSeq2Seq: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignKSW2(ref, query);
        }

        return result;
    }

    cigar::Cigar_t globalAlignSeq2Profile(const ProfileMatrix& ref,
                                const std::string& ref_string,
                              const std::string& query,
                              const anchor::Anchors& anchors)
    {
        align::AlignConfig cfg;
        align::AlignConfig first_cfg;
        first_cfg.flag = KSW_EZ_GENERIC_SC;

        const std::size_t ref_len = ref_string.size();
        const std::size_t qry_len = query.size();

        // 链化锚点：获取最佳链
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return globalAlignPSW(ref, query, cfg);
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
                                  std::size_t qry_start, std::size_t qry_end, align::AlignConfig seg_cfg,
                                  bool reverse_for_align = false) {
            // 边界裁剪
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            //const std::string seg_ref = ref_string.substr(ref_start, ref_end - ref_start);
            std::size_t seg_ref_len = ref_end - ref_start;
            ProfileMatrix seg_ref;
            seg_ref.len = static_cast<int>(seg_ref_len);
            seg_ref.dim = ref.dim;
            seg_ref.depth = ref.depth;
            if (seg_ref_len > 0) {
                const std::size_t offset = ref_start * static_cast<std::size_t>(ref.dim);
                const std::size_t count = seg_ref_len * static_cast<std::size_t>(ref.dim);
                seg_ref.prof.assign(ref.prof.begin() + static_cast<std::ptrdiff_t>(offset),
                                    ref.prof.begin() + static_cast<std::ptrdiff_t>(offset + count));
            }

            std::string seg_qry = query.substr(qry_start, qry_end - qry_start);
            const std::size_t seg_qry_len = seg_qry.size();

            // 仅用于尾段的质量优化：反向输入后求解，再把 CIGAR 顺序回正。
            // 注意：ref/query 角色不变，因此只需要反转 CIGAR 单元顺序，不需要互换 I/D。
            if (reverse_for_align) {
                std::reverse(seg_qry.begin(), seg_qry.end());

                // ProfileMatrix 是按“列块(dim)”线性存储，反向时必须按列翻转，
                // 否则会破坏单列内部 A/C/G/T/N 计数布局。
                if (seg_ref_len > 1) {
                    const std::size_t dim = static_cast<std::size_t>(seg_ref.dim);
                    std::vector<uint32_t> reversed_prof(seg_ref.prof.size(), 0U);
                    for (std::size_t dst_col = 0; dst_col < seg_ref_len; ++dst_col) {
                        const std::size_t src_col = (seg_ref_len - 1U) - dst_col;
                        const std::size_t src_off = src_col * dim;
                        const std::size_t dst_off = dst_col * dim;
                        std::copy(seg_ref.prof.begin() + static_cast<std::ptrdiff_t>(src_off),
                                  seg_ref.prof.begin() + static_cast<std::ptrdiff_t>(src_off + dim),
                                  reversed_prof.begin() + static_cast<std::ptrdiff_t>(dst_off));
                    }
                    seg_ref.prof.swap(reversed_prof);
                }
            }

            cigar::Cigar_t seg_cigar = globalAlignPSW(seg_ref, seg_qry, seg_cfg);
            if (reverse_for_align) {
                std::reverse(seg_cigar.begin(), seg_cigar.end());
            }

            // 用 CIGAR 反推消耗长度

            const std::size_t c_ref = cigar::getRefLength(seg_cigar);
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar);

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
#ifdef _DEBUG
                spdlog::warn("globalAlignSeq2Profile(seg): segment cigar mismatch (expected ref:{}/qry:{}, got ref:{}/qry:{}); forcing robust fallback for this segment",
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
        // 尾段启用“反向比对 + CIGAR 回正”，仅影响该段求解过程，不改变最终输出方向。
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg, true);

        // 最终一致性检查
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignSeq2Seq: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignPSW(ref, query, cfg);
        }

        return result;
    }

    cigar::Cigar_t globalAlignSeq2ProfileParallel(const ProfileMatrix& ref,
                                               const std::string& ref_string,
                                               const std::string& query,
                                               const anchor::Anchors& anchors,
                                               int thread)
    {
        align::AlignConfig cfg;
        align::AlignConfig first_cfg;
        first_cfg.flag = KSW_EZ_GENERIC_SC;

        const std::size_t ref_len = ref_string.size();
        const std::size_t qry_len = query.size();

        // 链化锚点：获取最佳链
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return globalAlignPSW(ref, query, cfg);
        }

        // 按 query 坐标排序
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        struct SegmentTask {
            std::size_t ref_start = 0;
            std::size_t ref_end = 0;
            std::size_t qry_start = 0;
            std::size_t qry_end = 0;
            align::AlignConfig seg_cfg{};
            bool reverse_for_align = false;
        };

        std::vector<SegmentTask> tasks;
        tasks.reserve(chain_anchors.size() * 2 + 2);

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto push_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                std::size_t qry_start, std::size_t qry_end,
                                align::AlignConfig seg_cfg,
                                bool reverse_for_align = false) {
            // 边界裁剪：保持和串行版 append_segment 完全一致
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            SegmentTask task;
            task.ref_start = ref_start;
            task.ref_end = ref_end;
            task.qry_start = qry_start;
            task.qry_end = qry_end;
            task.seg_cfg = seg_cfg;
            task.reverse_for_align = reverse_for_align;

            tasks.push_back(std::move(task));

            // 串行版 append_segment 最终一定把游标推进到 ref_end/qry_end：
            // 1) 正常情况下 c_ref == seg_ref_len, c_qry == seg_qry_len
            // 2) fallback 情况下显式 ref_pos = ref_end, qry_pos = qry_end
            ref_pos = ref_end;
            qry_pos = qry_end;
        };

        // 左端：起点到第一个锚点
        {
            const auto& first = chain_anchors.front();
            push_segment(ref_pos, first.pos_ref, qry_pos, first.pos_qry, first_cfg);
        }

        // 逐锚点：处理 span 和 gap
        for (std::size_t i = 0; i < chain_anchors.size(); ++i) {
            const auto& a = chain_anchors[i];

            const std::size_t a_ref_start = static_cast<std::size_t>(a.pos_ref);
            const std::size_t a_qry_start = static_cast<std::size_t>(a.pos_qry);
            const std::size_t a_ref_end = a_ref_start + static_cast<std::size_t>(a.span);
            const std::size_t a_qry_end = a_qry_start + static_cast<std::size_t>(a.span);

            push_segment(ref_pos, a_ref_end, qry_pos, a_qry_end, cfg);

            if (i + 1 < chain_anchors.size()) {
                const auto& b = chain_anchors[i + 1];
                push_segment(ref_pos, b.pos_ref, qry_pos, b.pos_qry, cfg);
            }
        }

        // 右端：最后一个锚点到末尾
        // 保留串行版逻辑：尾段启用“反向比对 + CIGAR 回正”
        push_segment(ref_pos, ref_len, qry_pos, qry_len, cfg, true);

        std::vector<cigar::Cigar_t> task_cigars(tasks.size());

    #ifdef _OPENMP
        const int use_threads = thread > 0 ? thread : omp_get_max_threads();

    #pragma omp parallel for default(none) shared(tasks, task_cigars, ref, query) num_threads(use_threads) schedule(static)
    #endif
        for (int task_idx = 0; task_idx < static_cast<int>(tasks.size()); ++task_idx) {
            const SegmentTask& task = tasks[static_cast<std::size_t>(task_idx)];

            const std::size_t seg_ref_len = task.ref_end - task.ref_start;

            ProfileMatrix seg_ref;
            seg_ref.len = static_cast<int>(seg_ref_len);
            seg_ref.dim = ref.dim;
            seg_ref.depth = ref.depth;

            if (seg_ref_len > 0) {
                const std::size_t offset = task.ref_start * static_cast<std::size_t>(ref.dim);
                const std::size_t count = seg_ref_len * static_cast<std::size_t>(ref.dim);

                seg_ref.prof.assign(
                    ref.prof.begin() + static_cast<std::ptrdiff_t>(offset),
                    ref.prof.begin() + static_cast<std::ptrdiff_t>(offset + count)
                );
            }

            std::string seg_qry = query.substr(task.qry_start, task.qry_end - task.qry_start);
            const std::size_t seg_qry_len = seg_qry.size();

            // 保留串行版尾段反向比对逻辑
            if (task.reverse_for_align) {
                std::reverse(seg_qry.begin(), seg_qry.end());

                // ProfileMatrix 是按“列块(dim)”线性存储，反向时必须按列翻转
                if (seg_ref_len > 1) {
                    const std::size_t dim = static_cast<std::size_t>(seg_ref.dim);
                    std::vector<uint32_t> reversed_prof(seg_ref.prof.size(), 0U);

                    for (std::size_t dst_col = 0; dst_col < seg_ref_len; ++dst_col) {
                        const std::size_t src_col = (seg_ref_len - 1U) - dst_col;
                        const std::size_t src_off = src_col * dim;
                        const std::size_t dst_off = dst_col * dim;

                        std::copy(
                            seg_ref.prof.begin() + static_cast<std::ptrdiff_t>(src_off),
                            seg_ref.prof.begin() + static_cast<std::ptrdiff_t>(src_off + dim),
                            reversed_prof.begin() + static_cast<std::ptrdiff_t>(dst_off)
                        );
                    }

                    seg_ref.prof.swap(reversed_prof);
                }
            }

            cigar::Cigar_t seg_cigar = globalAlignPSW(seg_ref, seg_qry, task.seg_cfg);

            if (task.reverse_for_align) {
                std::reverse(seg_cigar.begin(), seg_cigar.end());
            }

            const std::size_t c_ref = cigar::getRefLength(seg_cigar);
            const std::size_t c_qry = cigar::getQueryLength(seg_cigar);

            if (c_ref != seg_ref_len || c_qry != seg_qry_len) {
    #ifdef _DEBUG
                spdlog::warn(
                    "globalAlignSeq2Profile(seg): segment cigar mismatch "
                    "(expected ref:{}/qry:{}, got ref:{}/qry:{}); "
                    "forcing robust fallback for this segment",
                    seg_ref_len, seg_qry_len, c_ref, c_qry
                );
    #endif

                // 兜底策略：Query 全 I、Ref 全 D
                cigar::Cigar_t forced_cigar;

                if (seg_qry_len > 0) {
                    forced_cigar.push_back(
                        cigar::cigarToInt('I', static_cast<uint32_t>(seg_qry_len))
                    );
                }

                if (seg_ref_len > 0) {
                    forced_cigar.push_back(
                        cigar::cigarToInt('D', static_cast<uint32_t>(seg_ref_len))
                    );
                }

                task_cigars[static_cast<std::size_t>(task_idx)] = std::move(forced_cigar);
            } else {
                task_cigars[static_cast<std::size_t>(task_idx)] = std::move(seg_cigar);
            }
        }

        // 按原始分段顺序合并，保证结果不受并行执行顺序影响
        cigar::Cigar_t result;
        result.reserve(chain_anchors.size() * 2 + 2);

        for (const auto& seg_cigar : task_cigars) {
            cigar::appendCigar(result, seg_cigar);
        }

        // 最终一致性检查
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);

        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error(
                "globalAlignSeq2Seq: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                total_ref, ref_len, total_qry, qry_len
            );
            return globalAlignPSW(ref, query, cfg);
        }

        return result;
    }

} // namespace align
