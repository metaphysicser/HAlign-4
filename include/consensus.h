#ifndef HALIGN4_CONSENSUS_H
#define HALIGN4_CONSENSUS_H

#include <cstddef>
#include "utils.h"
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <string>
#include <vector>


// ---------------------------------------------------------------------------
// TopKLongestSelector
// 说明（中文）：
// - 用途：在单次流式扫描中维护“最长的 K 条序列”（Top-K by length），用于选择用于构建共识的候选集合。
// - 设计目标：低内存占用（只保存 K 条候选），单次遍历 (streaming) 输入即可得结果，时间复杂度 O(N log K)，
//   适合 N 很大但 K 相对较小的场景。
// - 稳定性：当长度相同时，优先保留更早出现的序列（stable tie-break）。这通过 order 字段（递增计数器）实现。
// - 并发：本类不是线程安全的；若在多线程环境中需要并发插入，请在外部对其加锁或采用每线程局部 TopK 后合并的策略。
// - 内存注意：保存的 SeqRecord 可能包含较长的序列，若 K 很大或序列很长，请注意内存占用；可替代地只保存索引或文件偏移。
// ---------------------------------------------------------------------------
class TopKLongestSelector
{
public:
    explicit TopKLongestSelector(std::size_t k = 0);

    void reset(std::size_t k);

    // 传值：调用方可选择复制传入或 std::move 传入，内部会 move 到堆里
    // 说明：使用传值语义可以让调用者决定是否移动 SeqRecord，从而避免不必要的复制。
    void consider(seq_io::SeqRecord rec);

    std::size_t size() const;
    std::size_t capacity() const;
    bool empty() const;

    // 取出结果（按长度降序；同长度按输入出现顺序升序），并清空内部状态
    // 说明：返回的是 SeqRecord 的移动语义，避免复制大量字符串。
    std::vector<seq_io::SeqRecord> takeSortedDesc();

private:
    struct Item
    {
        std::size_t len{0};
        std::uint64_t order{0};   // 输入顺序（用于稳定 tie-break）
        seq_io::SeqRecord rec;
    };

    // “更差/更小”判定：长度更短更差；长度相同则 order 更大（更晚出现）更差
    static bool worseThan(const Item& a, const Item& b);

    // 候选是否比当前最差（堆顶）更好
    static bool betterThan(const Item& cand, const Item& worst);

    void siftUp(std::size_t idx);
    void siftDown(std::size_t idx);

private:
    std::size_t k_{0};
    std::uint64_t order_counter_{0};

    // 自实现 min-heap：heap_[0] 永远是“最差”的那条（最短/同长最晚）
    std::vector<Item> heap_;
};

// ---------------------------------------------------------------------------
// consensus 命名空间：共识序列生成相关的数据结构与函数声明
// 这里包含：SiteCount（每个位点的碱基计数）、ConsensusJson（用于序列计数导出）、
// 基于字符到索引的快速映射表、以及共识生成/输出的接口。
// ---------------------------------------------------------------------------
namespace consensus
{
    // SiteCount：记录某个位点在一列中不同类别的计数
    // 字段说明：
    // - a,c,g,t,u: 对应碱基计数（U 支持 RNA/U），
    // - n: 未知碱基（N）计数，dash: gap ('-' 或 '.') 的计数
    // 性能提示：当前使用 uint32_t，如果你的数据集极大（每列计数可能超过 2^32），
    // 可以把类型改为 uint64_t，但会增加内存使用。
    struct SiteCount
    {
        std::uint32_t a = 0;
        std::uint32_t c = 0;
        std::uint32_t g = 0;
        std::uint32_t t = 0;
        std::uint32_t u = 0;
        std::uint32_t n = 0;
        std::uint32_t dash = 0;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("A", a),
               cereal::make_nvp("C", c),
               cereal::make_nvp("G", g),
               cereal::make_nvp("T", t),
               cereal::make_nvp("U", u),
               cereal::make_nvp("N", n),
               cereal::make_nvp("-", dash));
        }
    };

    // ConsensusJson：用于序列化输出共识相关统计信息（例如写入 JSON）
    // 字段说明：
    // - num_seqs: 参与统计的序列总数
    // - aln_len: 对齐长度（每个位点的计数向量长度）
    // - counts: 每个位点的 SiteCount 向量，长度为 aln_len
    struct ConsensusJson
    {
        std::uint64_t num_seqs = 0;
        std::uint64_t aln_len = 0;
        std::vector<SiteCount> counts;

        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(cereal::make_nvp("num_seqs", num_seqs),
               cereal::make_nvp("aln_len", aln_len),
               cereal::make_nvp("counts", counts));
        }
    };

    // ConsensusResult：一次 MSA 统计同时产出的三类数据。
    // - gap_seq: 保留 gap-majority 列的共识序列，用于后续按 MSA 列坐标比对。
    // - seq: 删除 gap_seq 中 gap 后的共识序列，用于 sketch/相似度等纯序列逻辑。
    // - counts: 每列计数，可用于构建 ProfileMatrix。
    struct ConsensusResult
    {
        std::string gap_seq;
        std::string seq;
        ConsensusJson counts;
    };

    // -------------------- 字符到索引的映射表 --------------------
    // 目的：在统计过程中需要把字符快速映射为 0..6 的索引（便于数组索引和分支最小化），
    // 使用一个 256 大小的查表（ASCII/unsigned char 范围）以常量时间完成映射。
    // 映射规则（idx）： 0->A, 1->C, 2->G, 3->T, 4->U, 5->N, 6->gap('-' 或 '.')
    // 性能说明：使用 inline constexpr 初始化保证在编译期生成常量数组，避免运行期初始化开销。
    inline constexpr std::array<std::uint8_t, 256> k_base_map = []() consteval {
        std::array<std::uint8_t, 256> m{};

        for (std::size_t i = 0; i < m.size(); ++i) m[i] = 5; // default N

        auto set = [&](unsigned char ch, std::uint8_t idx) { m[ch] = idx; };

        set((unsigned char)'A', 0); set((unsigned char)'a', 0);
        set((unsigned char)'C', 1); set((unsigned char)'c', 1);
        set((unsigned char)'G', 2); set((unsigned char)'g', 2);
        set((unsigned char)'T', 3); set((unsigned char)'t', 3);
        set((unsigned char)'U', 4); set((unsigned char)'u', 4);
        set((unsigned char)'N', 5); set((unsigned char)'n', 5);

        // gap：'-' 和 '.' 都当作 gap
        set((unsigned char)'-', 6);
        set((unsigned char)'.', 6);

        return m;
    }();

    // 简单包装：把单个字符映射到索引
    inline std::uint8_t mapBase(char ch)
    {
        return k_base_map[(unsigned char)ch];
    }

    // -------------------- 共识选择策略 --------------------
    // pickConsensusChar：给定某个位点的统计计数，选择一个最终共识碱基
    // 策略说明（当前实现建议）：
    // - 只在 A/C/G/T/U 五者之间选择（不会返回 N 或 gap），以保证共识序列尽量是可解析的核苷酸序列；
    // - 在计数相等时使用固定优先级（例如 A > C > G > T > U）来确保可复现性；
    // - 可选扩展：如果你希望在低覆盖位点返回 N，请修改策略以在总计数低于阈值时返回 'N'。
    char pickConsensusChar(const SiteCount& sc);

    // pickConsensusCharWithGap：当 gap 是严格最多的类别时返回 '-'，
    // 否则沿用 pickConsensusChar 的非 gap 共识选择策略。
    char pickConsensusCharWithGap(const SiteCount& sc);

    // writeConsensusFasta：将最终共识序列写入 FASTA（仅写一条 >consensus）
    // 注意：函数实现应保证输出目录存在（调用前可使用 file_io::ensureParentDirExists）
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq);

    // writeCountsJson：把 ConsensusJson 使用 cereal 写为 JSON 文件
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj);

    // -------------------- 共识生成接口 --------------------
    // generateConsensusResult：给定已对齐 FASTA，统计每个位点并同时生成：
    // 1) 带 gap 的共识序列；2) 去 gap 的共识序列；3) 每列计数。
    // out_fasta 写出去 gap 共识序列，保持与历史 consensus.fasta 用途兼容。
    ConsensusResult generateConsensusResult(const FilePath& aligned_fasta,
                                            const FilePath& out_fasta,
                                            const FilePath& out_json,
                                            std::uint64_t seq_limit,
                                            int thread,
                                            size_t batch_size = 4096);

    // generateConsensusSequence：给定已对齐的 FASTA 文件，统计每个位点并生成共识序列。
    // 参数说明：
    // - aligned_fasta: 已对齐的 FASTA（每个序列长度应一致为 aln_len）
    // - out_fasta: 写出共识序列的 FASTA 文件路径
    // - out_json: 写出统计计数的 JSON 文件路径
    // - seq_limit: 若非 0，可限制处理的序列数量（用于调试/抽样）
    // - thread: 期望使用的线程数（传入后函数会在内部设置 OpenMP 线程数或用于线程池规模）
    // 返回值：去 gap 共识序列字符串（便于内存中进一步处理或测试断言）
    // 性能提示：实现应该支持按批读取、线程本地累加(避免频繁原子更新)、SoA 布局以便向量化，以及合并阶段的并行化。
    // 另外建议在实现中记录耗时以便基准分析。
    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread,
                                                       size_t batch_size = 4096);

} // namespace consensus

#endif // HALIGN4_CONSENSUS_H
