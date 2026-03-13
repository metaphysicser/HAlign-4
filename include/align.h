// align.h - HAlign-4 序列比对模块核心接口
// 包括：CIGAR 操作、比对算法（KSW2/WFA2/MM2）、参考序列比对器

#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include "ksw2.h"
#include <unordered_map>
#include <filesystem>
#include <string>
#include <vector>
#include <functional>
#include "config.hpp"
#include "consensus.h"
#include "preprocess.h"

// CIGAR 操作：编码、解析、序列投影
// - 压缩格式：uint32_t (高28位长度 + 低4位操作符)
// - 操作码：0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
namespace cigar
{
    using CigarUnit = uint32_t;  // 单个 CIGAR 操作
    using Cigar_t = std::vector<CigarUnit>;  // CIGAR 序列

    // 编码/解码
    CigarUnit cigarToInt(char operation, uint32_t len);
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len);

    // 检测与转换
    bool hasInsertion(const Cigar_t& cigar);
    std::string cigarToString(const Cigar_t& cigar);
    Cigar_t stringToCigar(const std::string& cigar_str);

    // 序列投影对齐
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar);
    void delQueryToRefByCigar(std::string& query, const Cigar_t& cigar);

    // CIGAR 操作
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add);
    std::size_t getRefLength(const Cigar_t& cigar);
    std::size_t getQueryLength(const Cigar_t& cigar);
}

// 序列比对：KSW2、WFA2、锚点分段（MM2）
namespace align {
    // ------------------------------------------------------------------
    // 类型别名：种子（Seed）与种子命中（Seed Hit）
    // ------------------------------------------------------------------
    // 说明：
    // - 种子（Seed）是序列中的短 k-mer 或 minimizer，用于快速定位潜在的同源区域
    // - 种子命中（Seed Hit）记录种子在参考序列和查询序列中的位置
    // - 当前实现使用 minimizer 作为种子策略（高效、内存友好）
    // ------------------------------------------------------------------
    using SeedHit = minimizer::MinimizerHit;   // 单个种子命中：(ref_pos, query_pos, hash)
    using SeedHits = std::vector<SeedHit>;     // 种子命中列表（用于锚点定位）
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;  // 种子策略：minimizer

    // ------------------------------------------------------------------
    // DNA 字符到索引的映射表（编译期常量）
    // ------------------------------------------------------------------
    // 说明：将 DNA 字符（A/C/G/T/N，大小写不敏感）映射到 0-4 的索引
    // - 'A'/'a' -> 0
    // - 'C'/'c' -> 1
    // - 'G'/'g' -> 2
    // - 'T'/'t' -> 3
    // - 'N'/'n' 或其他 -> 4（未知碱基）
    // 用途：KSW2 需要将序列编码为整数数组才能进行比对
    // ------------------------------------------------------------------
    static constexpr uint8_t ScoreChar2Idx[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 0-15
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 16-31
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 32-47 (空格等)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 48-63 (数字)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 64-79  (@,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 80-95  (P,Q,R,S,T,U,V,W,X,Y,Z,...)
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,  // 96-111 (`,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,  // 112-127(p,q,r,s,t,u,v,w,x,y,z,...)
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 128-143
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 144-159
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 160-175
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 176-191
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 192-207
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 208-223
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 224-239
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4   // 240-255
    };

    // ------------------------------------------------------------------
    // DNA5 替换矩阵（编译期常量，用于 KSW2 全局/延伸比对）
    // ------------------------------------------------------------------
    // 说明：
    // 1. 矩阵维度：5x5，索引 0-3 为 A/C/G/T，索引 4 为 N（未知碱基/通配符）
    // 2. 一维展开：mat[i*5+j] 对应二维 mat[i][j]
    // 3. 评分规则：
    //    - 精确匹配（A-A, C-C, G-G, T-T）：+5 分（奖励正确配对）
    //    - 错配（A-C, A-G 等）：-4 分（惩罚碱基不匹配）
    //    - 涉及 N（未知碱基）：0 分（既不奖励也不惩罚）
    // 4. 参数平衡：
    //    - match/mismatch 比例约 5:4，与常见 Blast 参数（+1/-1 或 +5/-4）一致
    //    - 配合 gap_open=6, gap_extend=2 时，indel 代价约为 6+2k
    //    - 单个错配代价 9（从 +5 降到 -4），单个 1bp gap 代价 8，略偏好 gap
    // 5. 使用要求：
    //    - **必须**配合 KSW_EZ_GENERIC_SC flag 使用（启用完整替换矩阵）
    //    - 不带该 flag 时，KSW2 只支持简单 match/mismatch 模式，会忽略此矩阵
    // ------------------------------------------------------------------
    static constexpr int8_t dna5_simd_mat[25] = {
        // A   C   G   T   N
        5, -4, -4, -4,  0,  // A (i=0)
       -4,  5, -4, -4,  0,  // C (i=1)
       -4, -4,  5, -4,  0,  // G (i=2)
       -4, -4, -4,  5,  0,  // T (i=3)
        0,  0,  0,  0,  0   // N (i=4)
 };

    // ==================================================================
    // KSW2AlignConfig：KSW2 比对算法的配置结构体
    // ==================================================================
    // 说明：
    // KSW2 是一个高性能的序列比对算法（SSE/AVX 加速），支持全局、局部和延伸比对
    // 本结构体封装了 KSW2 所需的所有参数
    //
    // ==================================================================
    struct KSW2AlignConfig {
        const int8_t* mat = dna5_simd_mat;  // 替换矩阵（flattened 5x5，A/C/G/T/N）
        int alphabet_size = 5;              // 字母表大小（DNA5=5）
        int gap_open = 6;                   // gap open 罚分（正数，内部按罚分处理）
        int gap_extend = 2;                 // gap extend 罚分（正数）
        int end_bonus = 0;                  // 末端奖励（用于延伸比对，可为 0）
        int zdrop = -1;                     // Z-drop（-1 表示使用库默认/不启用）
        int band_width = -1;                // 带宽（-1 表示全矩阵，不限带宽）
        int flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT; // KSW2 标志：启用全替换矩阵等
    };


    // ==================================================================
    // 函数：auto_band - 自动估算 KSW2 带宽参数
    // ==================================================================
    // 功能：
    // 根据序列长度和预期的 indel（插入/删除）比率，自动估算合适的带宽
    // 用于加速长序列的 KSW2 比对（限制 DP 矩阵的计算范围）
    //
    // 参数：
    // @param qlen - query 序列长度
    // @param tlen - target（参考）序列长度
    // @param indel_rate - 预期的 indel 比率（默认 0.1 = 10%）
    //                     - 0.05：高相似度序列（如同种基因组内）
    //                     - 0.10：中等相似度（如近缘物种）
    //                     - 0.20：低相似度（如远缘物种）
    // @param margin - 安全边距（默认 200）
    //                 - 额外的缓冲区，防止真实 indel 超出带宽导致错误
    //
    // 返回：
    // 建议的带宽值（整数）
    //
    // ==================================================================
    //------------------------------------------- 带宽估计
    inline int auto_band(int qlen, int tlen,
        double indel_rate = 0.1,
        int    margin = 200)
    {
        // qlen/tlen 差异过大时不适合 banded DP：返回 -1 表示“禁用带宽限制”
        if ((double)std::abs(qlen - tlen) / (double)std::max(qlen, tlen) > 0.5)
        {
            return -1;
        }

        // 带宽≈预期 indel 规模 + 安全边距（经验公式，不影响正确性，只影响速度）
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));

    }

    // ==================================================================
    // 序列比对算法接口
    // ==================================================================
    // 说明：
    // 以下函数提供统一的序列比对接口，均返回 CIGAR 格式的比对结果
    // 输入序列可包含 A/C/G/T/N（大小写不敏感），其他字符视为 N
    // ==================================================================

    // ------------------------------------------------------------------
    // 函数：globalAlignKSW2 - KSW2 全局比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 KSW2 算法执行全局序列比对（Needleman-Wunsch 模式）
    // 适用场景：两条序列完整比对，考虑两端的所有碱基
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    //
    // 配置：
    // - 使用 dna5_simd_mat（5x5 替换矩阵）
    // - gap_open=6, gap_extend=2
    // - 全矩阵（无带宽限制）
    //
    // 性能：
    // - 时间复杂度：O(M * N)，M 和 N 分别为两序列长度
    // - 空间复杂度：O(M * N)（DP 矩阵）
    // - 加速：SSE/AVX 指令集加速（比普通 DP 快 4-8 倍）
    //
    // 示例：
    // globalAlignKSW2("ACGT", "AGT") -> "1M1D2M"（query 缺失了 ref 的 C）
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query, align::KSW2AlignConfig cfg);

    // ------------------------------------------------------------------
    // 函数：extendAlignKSW2 - KSW2 延伸比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 KSW2 算法执行延伸比对（extension alignment）
    // 适用场景：从种子（seed）位置向两侧延伸，直到得分下降过多或到达序列末端
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    // @param zdrop - Z-drop 剪枝阈值（默认 200）
    //                - 当比对得分下降超过 zdrop 时停止延伸
    //                - 值越大，延伸越完整但计算越慢
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述延伸区域的比对
    //
    // 配置：
    // - 使用 dna5_simd_mat（5x5 替换矩阵）
    // - gap_open=6, gap_extend=2
    // - end_bonus=5（奖励到达序列末端）
    // - 自动估算带宽（auto_band）
    //
    // 性能：
    // - 通常比全局比对快（因为剪枝和带宽限制）
    // - 适合长序列的局部比对（如基因组拼接中的重叠区域检测）
    //
    // 使用建议：
    // - 如果序列高度相似，可降低 zdrop（如 100）提高速度
    // - 如果需要完整比对，应使用 globalAlignKSW2 而非 extendAlignKSW2
    // ------------------------------------------------------------------
    cigar::Cigar_t extendAlignKSW2(const std::string& ref, const std::string& query, int zdrop = 200);

    // ------------------------------------------------------------------
    // 函数：globalAlignWFA2 - WFA2 全局比对
    // ------------------------------------------------------------------
    // 功能：
    // 使用 WFA2（Wavefront Alignment）算法执行全局序列比对
    // WFA2 特点：对于高相似度序列（indel 少），速度显著快于传统 DP 算法
    //
    // 参数：
    // @param ref - 参考序列（字符串，A/C/G/T/N）
    // @param query - 查询序列（字符串，A/C/G/T/N）
    //
    // 返回：
    // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref
    //
    // 性能：
    // - 时间复杂度：O(s * N)，s 为编辑距离，N 为序列长度
    //   - 对于高相似度序列（s << N），WFA2 远快于 KSW2
    //   - 对于低相似度序列（s ≈ N），WFA2 可能比 KSW2 慢
    // - 空间复杂度：O(s^2)（波前矩阵）
    //
    // 使用建议：
    // - 推荐用于高相似度序列（如同种基因组、近缘物种）
    // - 对于低相似度或未知相似度的序列，建议先用 MinHash 估算相似度
    // - 如果估算相似度 > 90%，使用 WFA2；否则使用 KSW2
    //
    // 示例：
    // globalAlignWFA2("ACGT", "AGT") -> "1M1D2M"
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);

    // cigar::Cigar_t extendAlignWFA2(const std::string& ref, const std::string& query, int zdrop = 200);

    // ------------------------------------------------------------------
    // 类型别名：比对函数类型
    // ------------------------------------------------------------------
    // 用于 globalAlignMM2 等函数，可以传入不同的比对算法（KSW2 或 WFA2）
    // 函数签名：接受两个序列（ref, query），返回 CIGAR
    using AlignFunc = std::function<cigar::Cigar_t(const std::string&, const std::string&)>;

    // ------------------------------------------------------------------
    // 函数：globalAlignMM2
    // 功能：基于锚点（anchors）的分段全局比对
    // 参数：
    //   - ref: 参考序列
    //   - query: 查询序列
    //   - anchors: 锚点集合（用于链化和分段）
    //   - align_func: 可选的比对函数（默认使用 globalAlignKSW2）
    //                 可以传入 globalAlignWFA2 或其他符合签名的比对函数
    // 返回：完整的 CIGAR 序列
    // 说明：
    //   - 使用锚点将序列分解为多个小片段，分别比对后拼接
    //   - 通过 align_func 参数可以灵活选择底层比对算法
    //   - 如果锚点为空或无效，会退化为 align_func 的全局比对
    // ------------------------------------------------------------------
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors);

    // ==================================================================
    // RefAligner 类：高性能参考序列比对器（多序列比对 MSA 引擎）
    // ==================================================================
    // 功能概述：
    // RefAligner 是 HAlign-4 的核心组件，负责将大量 query 序列比对到参考序列集合
    // 并生成最终的多序列比对（MSA）结果
    // 使用示例：
    // ------------------------------------------------------------------
    // // 1. 创建 RefAligner
    // Options opt = parseCommandLine(argc, argv);
    // RefAligner aligner(opt, "refs.fasta");
    //
    // // 2. 比对 query 序列
    // aligner.alignQueryToRef("queries.fasta", 5120);
    //
    // // 3. 合并结果生成 MSA
    // aligner.mergeAlignedResults("consensus_aligned.fasta", "mafft --auto");
    // ------------------------------------------------------------------
    // ==================================================================
    class RefAligner
    {
        public:
        // ------------------------------------------------------------------
        // 构造函数1：直接传入参数初始化
        //
        // 接口语义（与当前实现一致，且是调用方需要理解的部分）：
        // - 会读取 ref_fasta_path 中的参考序列集合，并为每条参考序列建立 sketch/minimizer 索引；
        // - 会生成一个“中心/共识序列”作为后续合并坐标系锚点：
        //   * keep_first_length=true  ：中心序列取第一条参考序列（ref_sequences[0]）；
        //   * keep_first_length=false ：中心序列由外部 MSA + generateConsensusSequence 生成。
        //
        // 参数：
        //   - work_dir：工作目录（中间文件与输出文件落盘位置）
        //   - ref_fasta_path：参考序列 FASTA
        //   - kmer_size/window_size/sketch_size：MinHash/Minimizer 参数
        //   - noncanonical：是否启用反向互补（DNA 常用 true）
        //   - threads：OpenMP 线程数（以及外部 MSA 可能用到）
        //   - msa_cmd：外部 MSA 命令模板（用于共识与插入序列对齐）
        //   - keep_first_length/keep_all_length：最终 MSA 列裁剪策略（见成员注释）
        // ------------------------------------------------------------------
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                   int kmer_size = 21, int window_size = 10,
                   int sketch_size = 2000, bool noncanonical = true,
                   int threads = 1, std::string msa_cmd = "",
                   bool keep_length = false,
                   bool enable_wfa = false);

        // ------------------------------------------------------------------
        // 构造函数2：基于 Options 结构体初始化（推荐方式）
        // ------------------------------------------------------------------
        // 参数：
        //   @param opt - Options 结构体（包含所有命令行参数和配置）
        //                - 自动提取：threads, kmer_size, window_size, sketch_size,
        //                  noncanonical, msa_cmd, keep_first_length, keep_all_length
        //   @param ref_fasta_path - 参考序列 FASTA 文件路径
        //
        // 工作流程：
        //   1. 读取参考序列文件到 ref_sequences
        //   2. 为每个参考序列生成 MinHash sketch（用于快速相似度估计）
        //   3. 为每个参考序列生成 Minimizer 索引（用于快速定位同源区域）
        //   4. 如果 opt.consensus_string 非空，直接使用；否则生成共识序列
        //   -------------------------
        RefAligner(const Options& opt, const FilePath& ref_fasta_path);

        // ==================================================================
        // 公共方法：比对与合并
        // ==================================================================

        // ------------------------------------------------------------------
        // 方法：alignQueryToRef - 并行比对 query 序列到参考序列
        // ------------------------------------------------------------------
        // 功能：
        // 读取 query FASTA 文件，将每条序列比对到最相似的参考序列
        // 输出多个 SAM 文件（每个线程一个），记录比对结果
        //
        // 参数：
        //   @param qry_fasta_path - query 序列 FASTA 文件路径
        //   @param batch_size - 批量读取的序列数量（默认 5120）
        //                       - 值越大，吞吐越高，但内存占用越多
        //                       - 建议：普通机器 1024-5120，高性能服务器 10000+
        //
        // ------------------------------------------------------------------
        // 说明：
        // - threads <= 0 表示使用 OpenMP 运行时默认线程数（例如由 OMP_NUM_THREADS 控制）
        // - batch_size 用于控制"流式读取"的批次大小，越大吞吐越高但占用内存更多
        void alignQueryToRef(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // ------------------------------------------------------------------
        // 方法：mergeAlignedResults - 合并所有比对结果生成最终 MSA
        // ------------------------------------------------------------------
        // 功能：
        // 将多个线程产生的 SAM 文件合并为一个统一的多序列比对（MSA）FASTA 文件
        // 包括处理插入序列的二次比对和坐标系统一
        //
        // 参数：
        //   @param msa_cmd - 外部 MSA 工具命令模板（用于对齐插入序列）
        //                    - 示例："mafft --auto {input} > {output}"
        //                    - 支持的工具：MAFFT, Muscle, Clustal Omega
        //   @param batch_size - 批处理大小：控制每批并行处理的序列数量
        //                       - 更大的 batch_size 可提高吞吐量，但占用更多内存
        //                       - 默认值 1000 适合大多数场景
        //   @param thread - 并行线程数：用于 OpenMP 并行处理 batch 内的序列
        //                   - 默认值 4，建议设置为 CPU 核心数
        // ------------------------------------------------------------------
        void mergeAlignedResults(const FilePath output, std::size_t batch_size = 25600);

        // ------------------------------------------------------------------
        // globalAlign - 全局序列比对（统一接口）
        // ------------------------------------------------------------------
        // 功能：
        // 根据两个序列的相似度自动选择合适的比对算法执行全局比对
        // 当前实现：直接使用 WFA2 算法（未来可扩展为相似度自适应策略）
        //
        // 参数：
        //   @param ref - 参考序列（字符串，A/C/G/T/N）
        //   @param query - 查询序列（字符串，A/C/G/T/N）
        //   @param similarity - 两个序列的相似度（0.0 到 1.0）
        //                       - 用于未来选择最优比对算法
        //                       - 当前版本未使用该参数，但保留接口兼容性
        //   @param ref_minimizer - 参考序列的 minimizer 索引（可选）
        //                          - 用于种子定位和锚点比对（未来扩展）
        //                          - 当前版本未使用，但保留接口用于优化
        //   @param query_minimizer - 查询序列的 minimizer 索引（可选）
        //                            - 用于种子定位和锚点比对（未来扩展）
        //                            - 当前版本未使用，但保留接口用于优化
        //
        // 返回：
        // CIGAR 操作序列（cigar::Cigar_t），描述 query 如何比对到 ref

        cigar::Cigar_t globalAlign(const std::string& ref,
                                          const std::string& query,
                                          double similarity,
                                          const SeedHits* ref_minimizer = nullptr,
                                          const SeedHits* query_minimizer = nullptr) const;

        // ------------------------------------------------------------------
        // 辅助函数：removeRefGapColumns
        // 功能：根据 ref_gap_pos 删除"参考为 gap 的那些列"（原地修改）
        //
        // 使用场景：
        // - 输入序列 seq 通常是"已经过 MSA 对齐"的序列（含 gap）
        // - ref_gap_pos 记录了参考（第一条序列）每一列是否为 gap
        // - 本函数会删除所有 ref_gap_pos[i]==true 的列，保留其余列
        //
        // 参数：
        //   - seq: 输入/输出序列（已对齐，包含 gap '-'）【原地修改】
        //   - ref_gap_pos: 参考（第一条序列）每一列是否为 gap；true 表示该列应被删除
        // ------------------------------------------------------------------
        static void removeRefGapColumns(
            std::string& seq,
            const std::vector<bool>& ref_gap_pos);


        private:
        // ------------------------------------------------------------------
        // 函数：alignOneQueryToRef
        // 功能：对单个 query 执行比对并写入 SAM 文件
        // 参数：
        //   - q: 待比对的查询序列
        //   - out: 普通输出文件的 writer（无插入或二次比对后无插入的序列）
        //   - out_insertion: 插入序列输出文件的 writer（二次比对后仍有插入的序列）
        // 说明：
        //   - 该函数不修改共享 reference 数据结构（线程安全）
        //   - out 和 out_insertion 由调用者管理生命周期（每线程独立）
        // ------------------------------------------------------------------
        void alignOneQueryToRef(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;


        // 辅助函数：写入SAM记录（选择正确的参考名称和输出文件）
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        // ------------------------------------------------------------------
        // 辅助函数：mergeConsensusAndSamToFasta
        // 功能：把 consensus_seq 与一组 SAM 文件中的序列合并写成一个 FASTA。
        //
        // keep 的语义（与当前实现一致）：
        // - keep=false：直接把 SAM 中的 query 序列原样写入 FASTA（不按 CIGAR 投影）。
        // - keep=true ：按 SAM CIGAR 将 query 投影到参考坐标（当前实现是“删除 query 的插入 I”，
        //              使其长度更接近参考/共识），用于后续 MSA/合并。
        //
        // 参数：
        //   - sam_paths：输入 SAM 文件路径列表
        //   - fasta_path：输出 FASTA 文件路径
        //   - keep：是否按 CIGAR 投影（见上）
        //   - line_width：FASTA 每行宽度（默认 80）
        // 返回：写入的总序列数（包含开头写入的 consensus_seq，至少为 1）
        // ------------------------------------------------------------------
        std::size_t mergeConsensusAndSamToFasta(
            const std::vector<FilePath>& sam_paths,
            const FilePath& fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map,
            bool keep = false,
            std::size_t line_width = 80
            ) const;

        // ------------------------------------------------------------------
        // 辅助函数：convertSamToFastaRecord
        // 功能：将单个 SAM 记录转换为 FASTA 记录，并根据 CIGAR 调整序列长度
        // ------------------------------------------------------------------
        void convertSamToFastaRecord(
            const seq_io::SamRecord& sam_rec,
            seq_io::SeqRecord& fasta_rec,
            const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
            std::size_t estimated_final_length) const;

        // ------------------------------------------------------------------
        // 辅助函数：processInsertionSequences
        // 功能：读取插入 SAM，合并为 FASTA，执行可选 MSA
        // 返回：插入序列 FASTA 文件路径
        // ------------------------------------------------------------------
        FilePath processInsertionSequences(
            const FilePath& result_dir,
            const FilePath& aligned_insertion_fasta,
            std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map) const;

        // ------------------------------------------------------------------
        // 辅助函数：writeConsensusAndReferences
        // 功能：从对齐文件读取并写入共识及参考序列
        // 返回：写入的序列数
        // ------------------------------------------------------------------
        std::size_t writeConsensusAndReferences(
            seq_io::SeqWriter& final_writer,
            const FilePath& consensus_aligned_file,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // 辅助函数：writeInsertionSequences
        // 功能：从对齐的插入文件读取并写入序列（跳过第一条共识）
        // 返回：写入的序列数
        // ------------------------------------------------------------------
        std::size_t writeInsertionSequences(
            seq_io::SeqWriter& final_writer,
            const FilePath& aligned_insertion_fasta,
            std::size_t& expected_length,
            bool& length_initialized,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // 辅助函数：processSamFileBatch
        // 功能：读取 SAM 批次，并行转换为 FASTA，串行写入输出
        // ------------------------------------------------------------------
        void processSamFileBatch(
            seq_io::SamReader& sam_reader,
            const std::size_t batch_size,
            seq_io::SeqWriter& final_writer,
            const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
            std::size_t estimated_final_length,
            std::size_t& expected_length,
            bool& length_initialized,
            std::size_t& seq_count,
            ProgressBar& progress) const;

        // ------------------------------------------------------------------
        // 辅助函数：parseAlignedReferencesToCigar
        // 功能：读取 MSA 对齐后的参考序列文件，生成每个序列的 CIGAR（只包含 M 和 D）
        //
        // 重要变更（接口约定）：
        // 1) 不再通过返回值返回 map，而是通过参数输出（避免大对象返回/移动，调用端更明确）
        // 2) 新增 ref_gap_pos：标记"参考序列对齐后的每一列是否为 gap（'-'）"
        //    - 这里的"参考序列"指该对齐文件中的第一条序列（通常是 consensus 或中心序列）
        //    - ref_gap_pos[i] == true  表示第 i 列参考为 gap
        //    - ref_gap_pos[i] == false 表示第 i 列参考为碱基
        //
        // 说明：
        // - 对于每个后续序列（参考/插入序列），只根据其自身字符是否为 '-' 来编码：碱基 -> M，gap -> D
        // - 不做碱基一致性校验（例如与第一条序列的列一致性），保持当前逻辑等价
        // ------------------------------------------------------------------
        void parseAlignedReferencesToCigar(
            const FilePath& aligned_fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
            std::vector<bool>& out_ref_gap_pos) const;


        // ==================================================================
        // 私有成员变量
        // ==================================================================

        // ------------------------------------------------------------------
        // 工作目录与数据路径
        // ------------------------------------------------------------------
        FilePath work_dir;  // 工作目录（存放中间文件和最终结果）

        // ------------------------------------------------------------------
        // 参考序列与索引
        // ------------------------------------------------------------------
        seq_io::SeqRecords ref_sequences;   // 参考序列集合（从 ref_fasta_path 读取）
        mash::Sketches ref_sketch;          // 每个参考序列的 MinHash sketch（用于快速相似度估计）
        std::vector<SeedHits> ref_minimizers;  // 每个参考序列的 Minimizer 索引（用于快速定位同源区域）

        // ------------------------------------------------------------------
        // 共识序列
        // ------------------------------------------------------------------
        // 说明：
        // - 共识序列是所有参考序列的"代表序列"
        // - 如果 keep_first_length=false，通过 MSA 工具生成（如 MAFFT）
        // - 如果 keep_first_length=true，直接使用 ref_sequences[0]
        // - 在 mergeAlignedResults 中作为坐标系的锚点
        // ------------------------------------------------------------------
        seq_io::SeqRecord consensus_seq;  // 共识序列（作为成员变量存储）
        mash::Sketch consensus_sketch;    // 共识序列的 MinHash sketch（用于二次比对时的相似度计算）
                                          // 说明：在构造函数中初始化一次，避免每次 alignOneQueryToRef 都重复计算
        SeedHits consensus_minimizer;     // 共识序列的 Minimizer 索引（用于二次比对时的种子定位）
                                          // 说明：在构造函数中初始化一次，避免每次 alignOneQueryToRef 都重复计算

        // ------------------------------------------------------------------
        // MinHash 和 Minimizer 参数
        // ------------------------------------------------------------------
        int kmer_size = 21;      // k-mer 大小（MinHash 和 Minimizer 共用）
        int window_size = 10;    // Minimizer 窗口大小
        int sketch_size = 2000;  // MinHash sketch 大小（哈希数量）
        int random_seed = 42;    // 随机种子（用于 MinHash 哈希函数）

        // ------------------------------------------------------------------
        // 线程与外部工具配置
        // ------------------------------------------------------------------
        int threads = 1;            // OpenMP 线程数（<=0 时通常表示让运行时决定；具体逻辑在 .cpp）
        std::string msa_cmd;        // 外部 MSA 命令模板（用于共识生成与插入序列 MSA）

        bool keep_length = false; // true：裁剪“共识为 gap 的列”，保持中心序列原始长度
        bool enable_wfa = false;  // true：允许使用 WFA 路径；默认 false 保持现有行为

        // ------------------------------------------------------------------
        // MinHash 计算选项
        // ------------------------------------------------------------------
        // 说明：
        // - noncanonical：是否使用非规范化的 k-mer（不区分正负链）
        //   * true：考虑 k-mer 的反向互补（用于 DNA 序列）
        //   * false：只考虑原始 k-mer（用于蛋白质或单链 RNA）
        // ------------------------------------------------------------------
        bool noncanonical = true;   // minimizer/sketch 是否启用反向互补（DNA 常用 true）

        // ------------------------------------------------------------------
        // 输出文件路径
        // ------------------------------------------------------------------
        // 说明：
        // - outs_path：每个线程的普通 SAM 输出文件路径
        //   * 长度 = threads
        //   * outs_path[i] = {work_dir}/data/aligned_thread_{i}.sam
        // - outs_with_insertion_path：每个线程的插入序列 SAM 输出文件路径
        //   * 长度 = threads
        //   * outs_with_insertion_path[i] = {work_dir}/data/aligned_insertion_thread_{i}.sam
        // ------------------------------------------------------------------
        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;
    };

} // namespace align

#endif //HALIGN4_ALIGN_H
