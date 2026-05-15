// align.h - HAlign-4 序列比对模块核心接口
// 包括：CIGAR 操作、比对算法（KSW2/WFA2/MM2）、参考序列比对器

#ifndef HALIGN4_ALIGN_H
#define HALIGN4_ALIGN_H
#include "utils.h"
#include "mash.h"
#include "seed.h"
#include "ksw2.h"
#include "psw.h"
#include <unordered_map>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <stdexcept>
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

    // 查询与转换
    bool hasInsertion(const Cigar_t& cigar);
    std::string cigarToString(const Cigar_t& cigar);
    Cigar_t stringToCigar(const std::string& cigar_str);

    // 按 CIGAR 将 query 映射到参考坐标
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar);
    void delQueryToRefByCigar(std::string& query, const Cigar_t& cigar);

    // CIGAR 拼接与长度统计
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add);
    std::size_t getRefLength(const Cigar_t& cigar);
    std::size_t getQueryLength(const Cigar_t& cigar);
}

// 序列比对接口：KSW2 / WFA2 / 锚点分段（MM2）
namespace align {
    // 种子命中类型（当前统一使用 minimizer）
    using SeedHit = minimizer::MinimizerHit;   // (ref_pos, query_pos, hash)
    using SeedHits = std::vector<SeedHit>;
    static constexpr seed::SeedKind kSeedKind = seed::SeedKind::minimizer;
	constexpr std::size_t kSketchBitsetRefThreshold = 10;
	constexpr unsigned kProfileEqualWeightShift = 10;
	constexpr std::uint32_t kProfileEqualWeightScale = 1U << kProfileEqualWeightShift;
    struct ProfileMatrix{
        int len;
        int dim;
        int depth;              /* profile 总序列数 */
        std::vector<uint32_t> prof;   /* 每列 dim 个计数；前 m 个通常是 residue/base 计数 */

        ProfileMatrix();

        // 从单条序列构造 profile：每列仅一个碱基计数为 1，其余为 0。
        // 约定 A/C/G/T/N -> 0/1/2/3/4，非法字符按 N 处理，保持与项目 DNA5 语义一致。
        explicit ProfileMatrix(const std::string& seq);

        explicit ProfileMatrix(const consensus::ConsensusJson& cj);

        static int baseIndex(char ch);

        static bool isGap(char ch);

        static ProfileMatrix fromAlignedSequence(const std::string& seq);

        // 从已对齐序列构建 profile：列坐标保持 MSA 坐标；gap 不进入 DNA5 计数，
        // depth 仍记录参与统计的序列总数，使 gap-rich 列在 profile 中表现为低覆盖列。
        static ProfileMatrix fromAlignedSequences(const std::vector<std::string>& aligned_sequences);

        // 从 consensus 统计计数构建 profile，避免 RefAligner 重新扫描 MSA 文件。
        // U 按 T 合并；gap 只通过 depth 与列总碱基数的差值体现，不占用额外维度。
        static ProfileMatrix fromConsensusCounts(const consensus::ConsensusJson& cj);
    };

	using ProfileMatrixMap = std::unordered_map<std::string, ProfileMatrix>;

    ProfileMatrix normalizeProfileEqualWeight(const ProfileMatrix& profile);
    ProfileMatrix combineProfilesEqualWeight(const std::vector<const ProfileMatrix*>& profiles);
    std::string profileToGappedSequence(const ProfileMatrix& profile);
    seq_io::SeqRecord reverseComplementRecord(const seq_io::SeqRecord& rec);


    // DNA 字符映射到 0..4（A/C/G/T/N，大小写不敏感；其他字符按 N）
    extern const uint8_t ScoreChar2Idx[256];
    // TODO 增加用户自定义矩阵
    // DNA5 替换矩阵（A/C/G/T/N）：match=+5，mismatch=-4，涉及 N 为 0。
    // 该矩阵需配合 KSW_EZ_GENERIC_SC 使用。
 //    static constexpr int8_t dna5_simd_mat[25] = {
 //        // A   C   G   T   N
 //        5, -4, -4, -4,  0,  // A (i=0)
 //       -4,  5, -4, -4,  0,  // C (i=1)
 //       -4, -4,  5, -4,  0,  // G (i=2)
 //       -4, -4, -4,  5,  0,  // T (i=3)
 //        0,  0,  0,  0,  0   // N (i=4)
 // };
     // mafft 风格
	extern const int8_t dna5_simd_mat[25];

    // KSW2 参数配置（默认值与当前实现一致）
    struct AlignConfig {
        const int8_t* mat = dna5_simd_mat;  // 5x5 替换矩阵（扁平化）
        int alphabet_size = 5;              // DNA5
        int gap_open = 10;                   // gap open 罚分
        int gap_extend = 2;                 // gap extend 罚分
        int end_bonus = 0;                  // 末端奖励
        int zdrop = -1;                     // -1 表示默认/不启用
        int band_width = -1;                // -1 表示不限制带宽
        int flag = KSW_EZ_RIGHT;
    };


    // 自动估计 band 宽度：长度差异过大时返回 -1（禁用 band）
    //------------------------------------------- 带宽估计
    int auto_band(int qlen, int tlen, double indel_rate = 0.1, int margin = 200);

    // 对齐接口统一返回 CIGAR；输入允许 A/C/G/T/N（其他字符按 N 处理）

    // KSW2 全局比对（Needleman-Wunsch）
    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);

    cigar::Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query, align::AlignConfig cfg);

    // KSW2 延伸比对（支持 zdrop 剪枝）
    cigar::Cigar_t extendAlignKSW2(const std::string& ref, const std::string& query, int zdrop = 200);

    // WFA2 全局比对：高相似度序列通常更快
    cigar::Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);

    cigar::Cigar_t globalAlignPSW(const ProfileMatrix& ref, const std::string& query, align::AlignConfig cfg);

    // 可注入的二元比对函数签名（ref, query）-> CIGAR
    using AlignFunc = std::function<cigar::Cigar_t(const std::string&, const std::string&)>;

    // 基于锚点的分段全局比对；锚点无效时退化为普通全局比对
    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  align::AlignConfig cfg = align::AlignConfig{});

    cigar::Cigar_t globalAlignSeq2Seq(const std::string& ref,
                                      const std::string& query,
                                      const anchor::Anchors& anchors,
                                      align::AlignConfig cfg = align::AlignConfig{});

    cigar::Cigar_t globalAlignSeq2Profile(const ProfileMatrix& ref,
                                    const std::string& ref_string,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  align::AlignConfig cfg = align::AlignConfig{});

    cigar::Cigar_t globalAlignSeq2ProfileParallel(const ProfileMatrix& ref,
                                    const std::string& ref_string,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  int thread = 1,
                                  align::AlignConfig cfg = align::AlignConfig{});


    struct MergeOptions {
        std::size_t batch_size = 25600;
        bool keep_length = false;
        bool write_reference = true;
        FilePath insertion_tsv_path;
        InsertionMergeMode insertion_merge_mode = InsertionMergeMode::reference_guided;
    };

    // 参考序列比对器：批量比对 query，并合并生成最终 MSA
    class RefAligner
    {
        public:
        // 直接参数构造：读取参考、构建索引、准备共识序列
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                   int kmer_size = 19, int window_size = 19,
                   int sketch_size = 30000, int profile_ref_kmer_len = 10,
                   bool noncanonical = true,
                   int threads = 1, std::string msa_cmd = "",
                   bool enable_wfa = false,
                   const FilePath& ref_aligned_path = FilePath(),
                   std::array<int8_t, 25> score_matrix = DEFAULT_DNA5_SCORE_MATRIX,
                   int gap_open = 10,
                   int gap_extend = 2,
                   int profile_ref_min = 15,
                   int profile_ref_max = 40,
                   double profile_ref_min_similarity = 0.7,
                   bool detect_reverse_complement = false);

        // Options 构造（推荐）
        RefAligner(const Options& opt, const FilePath& ref_fasta_path);

        // 并行将 query 比对到参考并输出线程独立 SAM
        void alignSeq2Seq(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // 并行将 query 比对到 profile，并输出线程独立 SAM
        void alignSeq2Profile(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // 合并 SAM 与中间结果，输出最终 MSA
        void mergeAlignedResults(const FilePath output, const MergeOptions& options);
        void mergeAlignedResults(const FilePath output, std::size_t batch_size = 25600);

        // 全局比对统一入口（保留 similarity/minimizer 参数以兼容后续策略）
        cigar::Cigar_t Seq2SeqWithAnchor(const std::string& ref,
                                          const std::string& query,
                                          double similarity,
                                          const SeedHits* ref_minimizer = nullptr,
                                          const SeedHits* query_minimizer = nullptr) const;

        cigar::Cigar_t Seq2ProfileWithAnchor(const ProfileMatrix& ref,
                                    const std::string& ref_string,
                                   const std::string& query,
                                   double similarity,
                                   int thread = 0,
                                   const SeedHits* ref_minimizer = nullptr,
                                   const SeedHits* query_minimizer = nullptr) const;

        // 删除“参考为 gap”的列（原地修改）
        static void removeRefGapColumns(
            std::string& seq,
            const std::vector<bool>& ref_gap_pos);


        private:
        void loadReference(const FilePath& ref_fasta_path, const FilePath& ref_aligned_path);
        void rebuildNormalizedReferenceProfiles();
        void refreshNormalizedReferenceProfiles(const std::vector<std::string>& ref_ids);

        // 单条 query 比对并写入 SAM（每线程 writer 由调用方管理）
        void alignOneQueryToRef(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;

        void alignOneQueryToProfile(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion,
                               cigar::Cigar_t& out_cigar,
                               std::vector<std::string>& out_ref_ids,
                               seq_io::SeqRecord& out_profile_query,
                               int thread = 0) const;

        // 根据一个 chunk 的对齐结果增量更新 profile 计数（不改变 profile 形状）
        void updateProfilesFromChunk(
            const std::vector<seq_io::SeqRecord>& chunk,
            const std::vector<cigar::Cigar_t>& cigar_chunk,
            const std::vector<std::vector<std::string>>& ref_ids_chunk);

        std::vector<mash::SketchMatch> selectProfileMatches(const mash::Sketch& query_sketch) const;

        // 将单条 query 按 CIGAR 投影到 profile 列并累加碱基计数，成功返回 true
        static bool applyCigarToProfile(
            const std::string& query_seq,
            const cigar::Cigar_t& cigar,
            ProfileMatrix& target_profile);

        // 写入一条 SAM 记录
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        // 解析对齐参考 FASTA：输出每条序列 CIGAR（M/D）与参考 gap 列标记
        // 不做碱基一致性校验，保持现有逻辑。
        void parseAlignedReferencesToCigar(
            const FilePath& aligned_fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
            std::vector<bool>& out_ref_gap_pos) const;

        AlignConfig makeAlignConfig() const;

        // 私有成员

        // 工作目录
        FilePath work_dir;

        // 参考序列与索引
        seq_io::SeqRecordMap ref_sequences;   // 参考序列集合
        mash::SketchMap ref_sketch;          // 每条参考序列的 MinHash sketch
        mash::SketchBitsetIndex ref_sketch_bitset_index;
        ProfileMatrixMap ref_profile;    // 参考序列的碱基计数 profile（按列存储，便于向量化）
        ProfileMatrixMap normalized_ref_profile; // 等权归一化后的 profile cache

        // 共识序列与索引（构造时预计算，避免重复计算）
        seq_io::SeqRecord consensus_gap_seq;
        seq_io::SeqRecord consensus_seq;
        mash::Sketch consensus_sketch;
        SeedHits consensus_minimizer;
        SeedHits consensus_gap_minimizer;
        ProfileMatrix consensus_profile;

        // MinHash / minimizer 参数
        int kmer_size = 21;
        int window_size = 10;
        int sketch_size = 2000;
        // 仅用于 profile reference sketch 构建的 k-mer 大小，与 minimizer 的 kmer_size 解耦。
        int profile_ref_kmer_len = 10;
        int random_seed = 42;

        // 并行与外部工具配置
        int threads = 1;            // OpenMP 线程数（<=0 时由运行时决定）
        std::string msa_cmd;        // 外部 MSA 命令模板

        bool enable_wfa = false;  // true：允许使用 WFA 路径
        bool profile_alignment_mode = false; // true：本轮 SAM 以 profile/带 gap 共识坐标为参考
        std::array<int8_t, 25> score_matrix = DEFAULT_DNA5_SCORE_MATRIX;
        int gap_open = 10;
        int gap_extend = 2;
        int profile_ref_min = 15;
        int profile_ref_max = 40;
        double profile_ref_min_similarity = 0.7;
        bool detect_reverse_complement = false;

        // 是否考虑反向互补
        bool noncanonical = true;

        // 每线程输出路径（普通/含插入）
        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;
    };

} // namespace align

#endif //HALIGN4_ALIGN_H
