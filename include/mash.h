#ifndef HALIGN4_MASH_H
#define HALIGN4_MASH_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "hash.h"
#include "bloom_filter.hpp"
#include <unordered_map>

namespace mash
{
    // nt4_table: 将 DNA 碱基字符映射到 0/1/2/3，其他字符映射到 4
    // A/a -> 0, C/c -> 1, G/g -> 2, T/t/U/u -> 3, 其他 -> 4
    // ASCII 码值：A=65, C=67, G=71, T=84, U=85, a=97, c=99, g=103, t=116, u=117
    inline constexpr std::uint8_t nt4_table[256] = {
        // 0-15
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 16-31
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 32-47
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 48-63
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 64-79: A=65->0, C=67->1, G=71->2
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        // 80-95: T=84->3, U=85->3
        4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
        // 96-111: a=97->0, c=99->1, g=103->2
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        // 112-127: t=116->3, u=117->3
        4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
        // 128-143
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 144-159
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 160-175
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 176-191
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 192-207
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 208-223
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 224-239
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 240-255
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
    };

    // ------------------------------------------------------------
    // Sketch (sorted unique bottom-k hashes)
    // ------------------------------------------------------------
    struct Sketch
    {
        std::size_t k = 0;
        bool noncanonical = true;
        std::vector<hash_t> hashes;

        std::size_t size() const noexcept { return hashes.size(); }
        bool empty() const noexcept { return hashes.empty(); }
    };
    using Sketches = std::vector<Sketch>;
	using SketchMap = std::unordered_map<std::string, Sketch>;  // id -> sketch

    struct SketchMatch
    {
        std::string id;
        double similarity = 0.0;
        bool found = false;
    };

    class SketchBitsetIndex
    {
    public:
        void build(const SketchMap& sketches, bool drop_hashes_present_in_all_refs = true);

        bool empty() const noexcept { return ids_.empty(); }
        bool usable() const noexcept { return !ids_.empty() && words_per_sketch_ != 0; }
        std::size_t size() const noexcept { return ids_.size(); }
        std::size_t dictionarySize() const noexcept { return dictionary_.size(); }

        SketchMatch findBest(const Sketch& query) const;
        std::vector<SketchMatch> findTopK(const Sketch& query,
                                          std::size_t min_count,
                                          std::size_t max_count,
                                          double similarity_ratio) const;

    private:
        std::size_t k_ = 0;
        bool has_k_ = false;
        std::size_t words_per_sketch_ = 0;
        std::vector<hash_t> dictionary_;
        std::unordered_map<hash_t, std::uint32_t> hash_to_bit_;
        std::vector<std::string> ids_;
        std::vector<std::uint32_t> sketch_hash_counts_;
        std::vector<std::uint64_t> bitsets_;
    };

    std::vector<hash_t> commonHashesAcrossSketches(const SketchMap& sketches);
    std::size_t removeCommonHashesFromSketches(SketchMap& sketches);

    // ------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------
    // 生成 MinHash sketch：对序列做 k-mer hash，然后取 bottom-k（sketch_size）并 sort+unique。
    // 注意：mash 模块与 seed/minimizer 无关；这里的 w 仅为兼容旧接口保留，当前实现不使用。
    Sketch sketchFromSequence(const std::string& seq,
                             std::size_t k,
                             std::size_t sketch_size,
                             bool noncanonical = true,
                             int seed = 0);

    bloom_filter filterFromSketch(const Sketch& sk, double false_positive_rate = 0.0001, int seed = 0xA5A5A5A5);


    // ------------------------------------------------------------
    // Similarity & distance
    // ------------------------------------------------------------

    // - both empty => 1.0
    // - one empty  => 0.0
    double jaccard(const Sketch& a, const Sketch& b);

    double jaccard(const bloom_filter& a, const Sketch& b);

    // d = - 1/k * ln( 2j/(1+j) )
    // - if j<=0 => +inf
    // - if j>=1 => 0
    double mashDistanceFromJaccard(double j, std::size_t k);

    // ANI ~= (2j/(1+j))^(1/k)
    // - clamp to [0,1]
    double aniFromJaccard(double j, std::size_t k);

    // ANI ~= exp(-d)
    // - clamp to [0,1]
    double aniFromMashDistance(double d);


    // ------------------------------------------------------------
    // Low-level helpers (sorted unique)
    // ------------------------------------------------------------

    std::size_t intersectionSizeSortedUnique(const std::vector<hash_t>& a,
                                             const std::vector<hash_t>& b) noexcept;

} // namespace mash

#endif //HALIGN4_MASH_H

