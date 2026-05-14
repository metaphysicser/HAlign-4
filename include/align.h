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

    struct ProfileMatrix{
        int len;
        int dim;
        int depth;              /* profile 总序列数 */
        std::vector<uint32_t> prof;   /* 每列 dim 个计数；前 m 个通常是 residue/base 计数 */

        ProfileMatrix() : len(0), dim(5), depth(0), prof() {}

        // 从单条序列构造 profile：每列仅一个碱基计数为 1，其余为 0。
        // 约定 A/C/G/T/N -> 0/1/2/3/4，非法字符按 N 处理，保持与项目 DNA5 语义一致。
        explicit ProfileMatrix(const std::string& seq) : len(static_cast<int>(seq.size())), dim(5), depth(seq.empty() ? 0 : 1), prof(static_cast<std::size_t>(len) * 5, 0U) {
            for (int i = 0; i < len; ++i) {
                const char ch = seq[static_cast<std::size_t>(i)];
                const int idx = baseIndex(ch);
                prof[static_cast<std::size_t>(i) * 5 + static_cast<std::size_t>(idx)] = 1U;
            }
        }

        explicit ProfileMatrix(const consensus::ConsensusJson& cj)
            : ProfileMatrix(fromConsensusCounts(cj))
        {
        }

        static int baseIndex(char ch)
        {
            switch (ch) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            case 'U': case 'u': return 3; // RNA/U 按 T 处理
            case 'N': case 'n': return 4;
            default: return 4;
            }
        }

        static bool isGap(char ch)
        {
            return ch == '-' || ch == '.';
        }

        static ProfileMatrix fromAlignedSequence(const std::string& seq)
        {
            if (seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("ProfileMatrix::fromAlignedSequence: alignment is too long");
            }

            ProfileMatrix pm;
            pm.len = static_cast<int>(seq.size());
            pm.dim = 5;
            pm.depth = 1;
            pm.prof.assign(seq.size() * static_cast<std::size_t>(pm.dim), 0U);

            for (std::size_t i = 0; i < seq.size(); ++i) {
                const char ch = seq[i];
                if (isGap(ch)) {
                    continue;
                }
                const std::size_t idx = static_cast<std::size_t>(baseIndex(ch));
                ++pm.prof[i * static_cast<std::size_t>(pm.dim) + idx];
            }

            return pm;
        }

        // 从已对齐序列构建 profile：列坐标保持 MSA 坐标；gap 不进入 DNA5 计数，
        // depth 仍记录参与统计的序列总数，使 gap-rich 列在 profile 中表现为低覆盖列。
        static ProfileMatrix fromAlignedSequences(const std::vector<std::string>& aligned_sequences)
        {
            if (aligned_sequences.empty()) {
                return ProfileMatrix();
            }
            if (aligned_sequences.size() == 1) {
                return fromAlignedSequence(aligned_sequences.front());
            }

            const std::size_t aln_len = aligned_sequences.front().size();
            if (aln_len > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("ProfileMatrix::fromAlignedSequences: alignment is too long");
            }
            if (aligned_sequences.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("ProfileMatrix::fromAlignedSequences: profile depth is too large");
            }

            ProfileMatrix pm;
            pm.len = static_cast<int>(aln_len);
            pm.dim = 5;
            pm.depth = static_cast<int>(aligned_sequences.size());
            pm.prof.assign(aln_len * static_cast<std::size_t>(pm.dim), 0U);

            for (const std::string& seq : aligned_sequences) {
                if (seq.size() != aln_len) {
                    throw std::runtime_error("ProfileMatrix::fromAlignedSequences: alignment length mismatch");
                }

                for (std::size_t i = 0; i < aln_len; ++i) {
                    const char ch = seq[i];
                    if (isGap(ch)) {
                        continue;
                    }
                    const std::size_t idx = static_cast<std::size_t>(baseIndex(ch));
                    ++pm.prof[i * static_cast<std::size_t>(pm.dim) + idx];
                }
            }

            return pm;
        }

        // 从 consensus 统计计数构建 profile，避免 RefAligner 重新扫描 MSA 文件。
        // U 按 T 合并；gap 只通过 depth 与列总碱基数的差值体现，不占用额外维度。
        static ProfileMatrix fromConsensusCounts(const consensus::ConsensusJson& cj)
        {
            if (cj.aln_len > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("ProfileMatrix::fromConsensusCounts: alignment is too long");
            }
            if (cj.num_seqs > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("ProfileMatrix::fromConsensusCounts: profile depth is too large");
            }
            if (cj.counts.size() != static_cast<std::size_t>(cj.aln_len)) {
                throw std::runtime_error("ProfileMatrix::fromConsensusCounts: counts length mismatch");
            }

            ProfileMatrix pm;
            pm.len = static_cast<int>(cj.aln_len);
            pm.dim = 5;
            pm.depth = static_cast<int>(cj.num_seqs);
            pm.prof.assign(static_cast<std::size_t>(pm.len) * static_cast<std::size_t>(pm.dim), 0U);

            for (std::size_t i = 0; i < cj.counts.size(); ++i) {
                const consensus::SiteCount& sc = cj.counts[i];
                const std::size_t off = i * static_cast<std::size_t>(pm.dim);
                pm.prof[off + 0] = sc.a;
                pm.prof[off + 1] = sc.c;
                pm.prof[off + 2] = sc.g;
                const std::uint64_t t_total = static_cast<std::uint64_t>(sc.t) + static_cast<std::uint64_t>(sc.u);
                pm.prof[off + 3] = static_cast<std::uint32_t>(t_total);
                pm.prof[off + 4] = sc.n;
            }

            return pm;
        }
    };

	using ProfileMatrixMap = std::unordered_map<std::string, ProfileMatrix>;


    // DNA 字符映射到 0..4（A/C/G/T/N，大小写不敏感；其他字符按 N）
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
	static constexpr int8_t dna5_simd_mat[25] = {
		// A   C   G   T   N
		 4, -2,  1, -2,  0,  // A (i=0)
	    -2,  4, -2,  1,  0,  // C (i=1)
	     1, -2,  4, -2,  0,  // G (i=2)
	    -2,  1, -2,  4,  0,  // T (i=3)
		 0,  0,  0,  0,  0   // N (i=4)
 };

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
    inline int auto_band(int qlen, int tlen,
        double indel_rate = 0.1,
        int    margin = 200)
    {
        // 长度差异过大时不适合 banded DP
        if ((double)std::abs(qlen - tlen) / (double)std::max(qlen, tlen) > 0.5)
        {
            return -1;
        }
        // 经验公式：预期 indel 规模 + 安全边距
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));

    }

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


    // 参考序列比对器：批量比对 query，并合并生成最终 MSA
    class RefAligner
    {
        public:
        // 直接参数构造：读取参考、构建索引、准备共识序列
        RefAligner(const FilePath& work_dir, const FilePath& ref_fasta_path,
                   int kmer_size = 19, int window_size = 19,
                   int sketch_size = 30000, int sketch_kmer_size = 21,
                   bool noncanonical = true,
                   int threads = 1, std::string msa_cmd = "",
                   bool keep_length = false,
                   bool enable_wfa = false,
                   const FilePath& ref_aligned_path = FilePath(),
                   std::array<int8_t, 25> score_matrix = DEFAULT_DNA5_SCORE_MATRIX,
                   int gap_open = 10,
                   int gap_extend = 2);

        // Options 构造（推荐）
        RefAligner(const Options& opt, const FilePath& ref_fasta_path);

        // 并行将 query 比对到参考并输出线程独立 SAM
        void alignSeq2Seq(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // 并行将 query 比对到 profile，并输出线程独立 SAM
        void alignSeq2Profile(const FilePath& qry_fasta_path, std::size_t batch_size = 25600);

        // 合并 SAM 与中间结果，输出最终 MSA
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
        // 单条 query 比对并写入 SAM（每线程 writer 由调用方管理）
        void alignOneQueryToRef(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion) const;

        void alignOneQueryToProfile(const seq_io::SeqRecord& q,
                               seq_io::SeqWriter& out,
                               seq_io::SeqWriter& out_insertion,
                               cigar::Cigar_t& out_cigar,
                               std::string& out_ref_id,
                               int thread = 0) const;

        // 根据一个 chunk 的对齐结果增量更新 profile 计数（不改变 profile 形状）
        void updateProfilesFromChunk(
            const std::vector<seq_io::SeqRecord>& chunk,
            const std::vector<cigar::Cigar_t>& cigar_chunk,
            const std::vector<std::string>& ref_id_chunk);

        // 将单条 query 按 CIGAR 投影到 profile 列并累加碱基计数，成功返回 true
        static bool applyCigarToProfile(
            const std::string& query_seq,
            const cigar::Cigar_t& cigar,
            ProfileMatrix& target_profile);

        // 写入一条 SAM 记录
        void writeSamRecord(const seq_io::SeqRecord& q, const cigar::Cigar_t& cigar,
                           std::string_view ref_name, seq_io::SeqWriter& out) const;

        // 共识 + SAM 合并为 FASTA
        // keep=false：原样写 query；keep=true：按 CIGAR 投影（当前会去除 query 的 I）
        std::size_t mergeConsensusAndSamToFasta(
            const std::vector<FilePath>& sam_paths,
            const FilePath& fasta_path,
            std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map,
            bool keep = false,
            std::size_t line_width = 80
            ) const;

        // 单条 SAM 转 FASTA，并按 CIGAR 调整长度
        void convertSamToFastaRecord(
            const seq_io::SamRecord& sam_rec,
            seq_io::SeqRecord& fasta_rec,
            const std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map,
            std::size_t estimated_final_length) const;

        // 处理插入序列：SAM -> FASTA -> 可选外部 MSA
        FilePath processInsertionSequences(
            const FilePath& result_dir,
            const FilePath& aligned_insertion_fasta,
            std::unordered_map<std::string, cigar::Cigar_t>& ref_aligned_map) const;

        // 写入共识与参考序列
        std::size_t writeConsensusAndReferences(
            seq_io::SeqWriter& final_writer,
            const FilePath& consensus_aligned_file,
            ProgressBar& progress) const;

        // 写入插入序列（跳过首条共识）
        std::size_t writeInsertionSequences(
            seq_io::SeqWriter& final_writer,
            const FilePath& aligned_insertion_fasta,
            std::size_t& expected_length,
            bool& length_initialized,
            ProgressBar& progress) const;

        // 批量读取 SAM，并行转换 FASTA，串行写出
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
        // 仅用于 sketch 构建的 k-mer 大小，与 minimizer 的 kmer_size 解耦。
        int sketch_kmer_size = 21;
        int random_seed = 42;

        // 并行与外部工具配置
        int threads = 1;            // OpenMP 线程数（<=0 时由运行时决定）
        std::string msa_cmd;        // 外部 MSA 命令模板

        bool keep_length = false; // true：裁剪“共识为 gap”的列
        bool enable_wfa = false;  // true：允许使用 WFA 路径
        bool profile_alignment_mode = false; // true：本轮 SAM 以 profile/带 gap 共识坐标为参考
        std::array<int8_t, 25> score_matrix = DEFAULT_DNA5_SCORE_MATRIX;
        int gap_open = 10;
        int gap_extend = 2;

        // 是否考虑反向互补
        bool noncanonical = true;

        // 每线程输出路径（普通/含插入）
        std::vector<FilePath> outs_path;
        std::vector<FilePath> outs_with_insertion_path;
    };

} // namespace align

#endif //HALIGN4_ALIGN_H
