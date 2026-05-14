#include "mash.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace mash {
namespace detail {
    std::size_t bitsetAndPopcountSse2(const std::uint64_t* a,
                                      const std::uint64_t* b,
                                      std::size_t words) noexcept;

    std::size_t bitsetAndPopcountScalar(const std::uint64_t* a,
                                        const std::uint64_t* b,
                                        std::size_t words) noexcept
    {
        std::size_t count = 0;
        for (std::size_t i = 0; i < words; ++i) {
            count += static_cast<std::size_t>(std::popcount(a[i] & b[i]));
        }
        return count;
    }

    std::size_t bitsetAndPopcount(const std::uint64_t* a,
                                  const std::uint64_t* b,
                                  std::size_t words) noexcept
    {
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
        return bitsetAndPopcountSse2(a, b, words);
#else
        return bitsetAndPopcountScalar(a, b, words);
#endif
    }
} // namespace detail

    std::vector<hash_t> commonHashesAcrossSketches(const SketchMap& sketches)
    {
        std::vector<hash_t> common;
        const std::size_t ref_count = sketches.size();
        if (ref_count <= 1) {
            return common;
        }

        std::size_t total_hashes = 0;
        for (const auto& [id, sketch] : sketches) {
            (void)id;
            total_hashes += sketch.hashes.size();
        }

        std::unordered_map<hash_t, std::uint32_t> counts;
        counts.reserve(total_hashes);
        for (const auto& [id, sketch] : sketches) {
            (void)id;
            for (const hash_t h : sketch.hashes) {
                ++counts[h];
            }
        }

        common.reserve(counts.size());
        for (const auto& [h, count] : counts) {
            if (static_cast<std::size_t>(count) == ref_count) {
                common.push_back(h);
            }
        }
        std::sort(common.begin(), common.end());
        return common;
    }

    std::size_t removeCommonHashesFromSketches(SketchMap& sketches)
    {
        const std::vector<hash_t> common = commonHashesAcrossSketches(sketches);
        if (common.empty()) {
            return 0;
        }

        for (auto& [id, sketch] : sketches) {
            (void)id;
            std::vector<hash_t> filtered;
            filtered.reserve(sketch.hashes.size());
            std::set_difference(sketch.hashes.begin(), sketch.hashes.end(),
                                common.begin(), common.end(),
                                std::back_inserter(filtered));
            sketch.hashes.swap(filtered);
        }

        return common.size();
    }

    void SketchBitsetIndex::build(const SketchMap& sketches,
                                  bool drop_hashes_present_in_all_refs)
    {
        k_ = 0;
        has_k_ = false;
        words_per_sketch_ = 0;
        dictionary_.clear();
        hash_to_bit_.clear();
        ids_.clear();
        sketch_hash_counts_.clear();
        bitsets_.clear();

        if (sketches.empty()) {
            return;
        }

        ids_.reserve(sketches.size());
        for (const auto& [id, sketch] : sketches) {
            if (!has_k_) {
                k_ = sketch.k;
                has_k_ = true;
            } else if (sketch.k != k_) {
                throw std::invalid_argument("SketchBitsetIndex::build: mismatched sketch k");
            }
            ids_.push_back(id);
        }
        std::sort(ids_.begin(), ids_.end());

        std::size_t total_hashes = 0;
        for (const auto& [id, sketch] : sketches) {
            (void)id;
            total_hashes += sketch.hashes.size();
        }

        std::unordered_map<hash_t, std::uint32_t> counts;
        counts.reserve(total_hashes);
        for (const auto& [id, sketch] : sketches) {
            (void)id;
            for (const hash_t h : sketch.hashes) {
                ++counts[h];
            }
        }

        dictionary_.reserve(counts.size());
        const std::uint32_t ref_count = static_cast<std::uint32_t>(sketches.size());
        for (const auto& [h, count] : counts) {
            if (drop_hashes_present_in_all_refs && count == ref_count) {
                continue;
            }
            dictionary_.push_back(h);
        }
        std::sort(dictionary_.begin(), dictionary_.end());

        if (dictionary_.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::runtime_error("SketchBitsetIndex::build: dictionary is too large");
        }

        hash_to_bit_.reserve(dictionary_.size());
        for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(dictionary_.size()); ++i) {
            hash_to_bit_.emplace(dictionary_[i], i);
        }

        words_per_sketch_ = (dictionary_.size() + 63U) / 64U;
        if (words_per_sketch_ == 0) {
            return;
        }

        sketch_hash_counts_.assign(ids_.size(), 0);
        bitsets_.assign(ids_.size() * words_per_sketch_, 0);

        for (std::size_t row = 0; row < ids_.size(); ++row) {
            const auto sketch_it = sketches.find(ids_[row]);
            if (sketch_it == sketches.end()) {
                throw std::runtime_error("SketchBitsetIndex::build: id disappeared while building index");
            }

            std::uint64_t* bits = bitsets_.data() + row * words_per_sketch_;
            std::uint32_t retained = 0;
            for (const hash_t h : sketch_it->second.hashes) {
                const auto bit_it = hash_to_bit_.find(h);
                if (bit_it == hash_to_bit_.end()) {
                    continue;
                }
                const std::uint32_t bit = bit_it->second;
                bits[bit >> 6U] |= (1ULL << (bit & 63U));
                ++retained;
            }
            sketch_hash_counts_[row] = retained;
        }
    }

    SketchMatch SketchBitsetIndex::findBest(const Sketch& query) const
    {
        SketchMatch best;
        if (ids_.empty()) {
            return best;
        }
        if (has_k_ && query.k != k_) {
            throw std::invalid_argument("SketchBitsetIndex::findBest: mismatched sketch k");
        }
        if (words_per_sketch_ == 0 || dictionary_.empty()) {
            best.id = ids_.front();
            best.found = true;
            return best;
        }

        thread_local std::vector<std::uint64_t> query_words;
        query_words.assign(words_per_sketch_, 0);

        std::uint32_t query_hash_count = 0;
        for (const hash_t h : query.hashes) {
            const auto bit_it = hash_to_bit_.find(h);
            if (bit_it == hash_to_bit_.end()) {
                continue;
            }
            const std::uint32_t bit = bit_it->second;
            query_words[bit >> 6U] |= (1ULL << (bit & 63U));
            ++query_hash_count;
        }

        double best_score = -1.0;
        for (std::size_t row = 0; row < ids_.size(); ++row) {
            const std::uint64_t* ref_bits = bitsets_.data() + row * words_per_sketch_;
            const std::size_t inter = detail::bitsetAndPopcount(
                query_words.data(), ref_bits, words_per_sketch_);
            const std::uint32_t denom = std::min(query_hash_count, sketch_hash_counts_[row]);
            const double score = denom == 0 ? 0.0 :
                static_cast<double>(inter) / static_cast<double>(denom);

            if (!best.found || score > best_score) {
                best.id = ids_[row];
                best.similarity = score;
                best.found = true;
                best_score = score;
            }
        }

        return best;
    }
} // namespace mash
