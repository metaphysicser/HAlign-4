#include <bit>
#include <cstdint>
#include <cstddef>

#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#include <emmintrin.h>
#endif

namespace mash::detail {
    std::size_t bitsetAndPopcountScalar(const std::uint64_t* a,
                                        const std::uint64_t* b,
                                        std::size_t words) noexcept;

    std::size_t bitsetAndPopcountSse2(const std::uint64_t* a,
                                      const std::uint64_t* b,
                                      std::size_t words) noexcept
    {
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
        std::size_t count = 0;
        std::size_t i = 0;
        alignas(16) std::uint64_t tmp[2];

        for (; i + 1 < words; i += 2) {
            const __m128i av = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
            const __m128i bv = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i));
            const __m128i both = _mm_and_si128(av, bv);
            _mm_store_si128(reinterpret_cast<__m128i*>(tmp), both);
            count += static_cast<std::size_t>(std::popcount(tmp[0]));
            count += static_cast<std::size_t>(std::popcount(tmp[1]));
        }

        if (i < words) {
            count += static_cast<std::size_t>(std::popcount(a[i] & b[i]));
        }
        return count;
#else
        return bitsetAndPopcountScalar(a, b, words);
#endif
    }
} // namespace mash::detail
