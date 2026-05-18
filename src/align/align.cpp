#include "align.h"
#include "ksw2.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>

extern "C" {
#include "alignment/cigar.h"
#include "wavefront/wavefront_align.h"
}

// 序列比对算法封装：KSW2 / WFA2
// - 统一返回 cigar::Cigar_t（压缩格式）
// - 支持全局比对、延伸比对、锚点分段比对

namespace align
{
    namespace {

        static constexpr std::uint64_t kProfileBaseThresholdNum = 7;
        static constexpr std::uint64_t kProfileBaseThresholdDen = 10;
        static constexpr std::uint32_t kMaxCigarOpLen = (1U << 28) - 1U;

        struct MatchBlock
        {
            std::size_t ref_start = 0;
            std::size_t qry_start = 0;
            std::size_t len = 0;
        };

        char complementBase(char ch)
        {
            switch (ch) {
                case 'A': case 'a': return 'T';
                case 'C': case 'c': return 'G';
                case 'G': case 'g': return 'C';
                case 'T': case 't': return 'A';
                case 'U': case 'u': return 'A';
                case '-': case '.': return ch;
                default: return 'N';
            }
        }

        char canonicalDnaBase(char ch)
        {
            switch (ch) {
                case 'A': case 'a': return 'A';
                case 'C': case 'c': return 'C';
                case 'G': case 'g': return 'G';
                case 'T': case 't':
                case 'U': case 'u': return 'T';
                default: return '\0';
            }
        }

        bool simpleBasesMatch(char ref_base, char query_base)
        {
            const char qry = canonicalDnaBase(query_base);
            return ref_base != '\0' && qry != '\0' && ref_base == qry;
        }

        char profileConsensusBaseAt(const ProfileMatrix& profile, std::size_t pos)
        {
            if (profile.len < 0 || profile.dim < 5 || profile.depth <= 0 ||
                pos >= static_cast<std::size_t>(profile.len)) {
                return '\0';
            }

            const std::size_t dim = static_cast<std::size_t>(profile.dim);
            const std::size_t off = pos * dim;
            if (off + 4 >= profile.prof.size()) {
                return '\0';
            }

            static constexpr char bases[4] = {'A', 'C', 'G', 'T'};
            std::uint32_t best_count = 0;
            std::size_t best_idx = 0;
            for (std::size_t idx = 0; idx < 4; ++idx) {
                const std::uint32_t count = profile.prof[off + idx];
                if (count > best_count) {
                    best_count = count;
                    best_idx = idx;
                }
            }

            const std::uint64_t lhs =
                static_cast<std::uint64_t>(best_count) * kProfileBaseThresholdDen;
            const std::uint64_t rhs =
                static_cast<std::uint64_t>(profile.depth) * kProfileBaseThresholdNum;
            return lhs > rhs ? bases[best_idx] : '\0';
        }

        void appendCigarOp(cigar::Cigar_t& result, char op, std::size_t len)
        {
            while (len > 0) {
                if (!result.empty()) {
                    char last_op = '\0';
                    std::uint32_t last_len = 0;
                    cigar::intToCigar(result.back(), last_op, last_len);
                    if (last_op == op && last_len < kMaxCigarOpLen) {
                        const std::size_t add =
                            std::min<std::size_t>(len, kMaxCigarOpLen - last_len);
                        result.back() = cigar::cigarToInt(
                            op, last_len + static_cast<std::uint32_t>(add));
                        len -= add;
                        continue;
                    }
                }

                const std::size_t chunk = std::min<std::size_t>(len, kMaxCigarOpLen);
                result.push_back(cigar::cigarToInt(op, static_cast<std::uint32_t>(chunk)));
                len -= chunk;
            }
        }

        template <typename RefBaseAt>
        std::vector<MatchBlock> buildExtendedMatchBlocks(const anchor::Anchors& chain_anchors,
                                                         std::size_t ref_len,
                                                         std::size_t qry_len,
                                                         const std::string& query,
                                                         RefBaseAt ref_base_at)
        {
            std::vector<MatchBlock> blocks;
            blocks.reserve(chain_anchors.size());

            std::size_t prev_ref_end = 0;
            std::size_t prev_qry_end = 0;

            for (const anchor::Anchor& anchor : chain_anchors) {
                std::size_t ref_start = std::min<std::size_t>(anchor.pos_ref, ref_len);
                std::size_t qry_start = std::min<std::size_t>(anchor.pos_qry, qry_len);

                std::size_t shift = 0;
                if (ref_start < prev_ref_end) {
                    shift = std::max(shift, prev_ref_end - ref_start);
                }
                if (qry_start < prev_qry_end) {
                    shift = std::max(shift, prev_qry_end - qry_start);
                }

                if (shift > ref_len - ref_start || shift > qry_len - qry_start) {
                    continue;
                }
                ref_start += shift;
                qry_start += shift;

                std::size_t max_len = static_cast<std::size_t>(anchor.span);
                max_len = std::min(max_len, ref_len - ref_start);
                max_len = std::min(max_len, qry_len - qry_start);
                if (shift >= static_cast<std::size_t>(anchor.span) || max_len == 0) {
                    continue;
                }

                std::size_t len = 0;
                while (len < max_len &&
                       simpleBasesMatch(ref_base_at(ref_start + len), query[qry_start + len])) {
                    ++len;
                }
                if (len == 0) {
                    continue;
                }

                while (ref_start > prev_ref_end && qry_start > prev_qry_end &&
                       simpleBasesMatch(ref_base_at(ref_start - 1), query[qry_start - 1])) {
                    --ref_start;
                    --qry_start;
                    ++len;
                }

                while (ref_start + len < ref_len &&
                       qry_start + len < qry_len &&
                       simpleBasesMatch(ref_base_at(ref_start + len), query[qry_start + len])) {
                    ++len;
                }

                blocks.push_back(MatchBlock{ref_start, qry_start, len});
                prev_ref_end = ref_start + len;
                prev_qry_end = qry_start + len;
            }

            return blocks;
        }

        ProfileMatrix sliceProfile(const ProfileMatrix& profile,
                                   std::size_t ref_start,
                                   std::size_t ref_end)
        {
            if (profile.dim <= 0 || profile.len < 0) {
                throw std::runtime_error("sliceProfile: invalid profile shape");
            }

            const std::size_t profile_len = static_cast<std::size_t>(profile.len);
            ref_start = std::min(ref_start, profile_len);
            ref_end = std::min(ref_end, profile_len);
            if (ref_end < ref_start) {
                ref_end = ref_start;
            }

            ProfileMatrix seg_ref;
            const std::size_t seg_ref_len = ref_end - ref_start;
            seg_ref.len = static_cast<int>(seg_ref_len);
            seg_ref.dim = profile.dim;
            seg_ref.depth = profile.depth;

            if (seg_ref_len > 0) {
                const std::size_t dim = static_cast<std::size_t>(profile.dim);
                const std::size_t offset = ref_start * dim;
                const std::size_t count = seg_ref_len * dim;
                if (offset + count > profile.prof.size()) {
                    throw std::runtime_error("sliceProfile: profile data is truncated");
                }
                seg_ref.prof.assign(
                    profile.prof.begin() + static_cast<std::ptrdiff_t>(offset),
                    profile.prof.begin() + static_cast<std::ptrdiff_t>(offset + count));
            }

            return seg_ref;
        }

        std::string kswCompatibleReference(std::string ref)
        {
            for (char& ch : ref) {
                if (ProfileMatrix::isGap(ch) || canonicalDnaBase(ch) == '\0') {
                    ch = 'N';
                }
            }
            return ref;
        }

        cigar::Cigar_t extendAlignKSW2Configured(const std::string& ref,
                                                 const std::string& query,
                                                 align::AlignConfig cfg)
        {
            if (ref.empty() || query.empty()) {
                return globalAlignKSW2(ref, query, cfg);
            }

            std::vector<std::uint8_t> ref_enc(ref.size());
            std::vector<std::uint8_t> qry_enc(query.size());
            for (std::size_t i = 0; i < ref.size(); ++i) {
                ref_enc[i] = align::ScoreChar2Idx[static_cast<std::uint8_t>(ref[i])];
            }
            for (std::size_t i = 0; i < query.size(); ++i) {
                qry_enc[i] = align::ScoreChar2Idx[static_cast<std::uint8_t>(query[i])];
            }

            cfg.zdrop = cfg.zdrop >= 0 ? cfg.zdrop : 200;
            cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
            cfg.end_bonus = cfg.end_bonus != 0 ? cfg.end_bonus : 50;
            if (cfg.band_width == AUTO_BAND_WIDTH) {
                cfg.band_width = align::auto_band(
                    static_cast<int>(ref.size()),
                    static_cast<int>(query.size()));
            }

            ksw_extz_t ez{};
            ksw_extz2_sse(nullptr,
                static_cast<int>(qry_enc.size()), qry_enc.data(),
                static_cast<int>(ref_enc.size()), ref_enc.data(),
                cfg.alphabet_size, cfg.mat,
                cfg.gap_open, cfg.gap_extend,
                cfg.band_width, cfg.zdrop, cfg.end_bonus,
                cfg.flag, &ez);

            cigar::Cigar_t out;
            out.reserve(ez.n_cigar);
            for (int i = 0; i < ez.n_cigar; ++i) {
                out.push_back(ez.cigar[i]);
            }

            free(ez.cigar);
            return out;
        }

        cigar::Cigar_t alignTerminalKSW2(const std::string& ref_segment,
                                         const std::string& query_segment,
                                         align::AlignConfig cfg,
                                         bool reverse_from_anchor)
        {
            std::string ref = kswCompatibleReference(ref_segment);
            std::string query = query_segment;

            if (ref.empty() || query.empty()) {
                return globalAlignKSW2(ref, query, cfg);
            }

            if (!reverse_from_anchor) {
                cigar::Cigar_t ext = extendAlignKSW2Configured(ref, query, cfg);
                const std::size_t c_ref = cigar::getRefLength(ext);
                const std::size_t c_qry = cigar::getQueryLength(ext);
                if (c_ref <= ref.size() && c_qry <= query.size() && !ext.empty()) {
                    if (c_ref < ref.size() || c_qry < query.size()) {
                        cigar::appendCigar(ext, globalAlignKSW2(
                            ref.substr(c_ref), query.substr(c_qry), cfg));
                    }
                    return ext;
                }
                return globalAlignKSW2(ref, query, cfg);
            }

            std::string rev_ref = ref;
            std::string rev_query = query;
            std::reverse(rev_ref.begin(), rev_ref.end());
            std::reverse(rev_query.begin(), rev_query.end());

            cigar::Cigar_t ext = extendAlignKSW2Configured(rev_ref, rev_query, cfg);
            const std::size_t c_ref = cigar::getRefLength(ext);
            const std::size_t c_qry = cigar::getQueryLength(ext);
            if (c_ref > ref.size() || c_qry > query.size() || ext.empty()) {
                return globalAlignKSW2(ref, query, cfg);
            }

            const std::size_t prefix_ref_len = ref.size() - c_ref;
            const std::size_t prefix_qry_len = query.size() - c_qry;
            cigar::Cigar_t out = globalAlignKSW2(
                ref.substr(0, prefix_ref_len),
                query.substr(0, prefix_qry_len),
                cfg);
            std::reverse(ext.begin(), ext.end());
            cigar::appendCigar(out, ext);
            return out;
        }
    } // namespace

    const uint8_t ScoreChar2Idx[256] = {
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

    const int8_t dna5_simd_mat[25] = {
        // A   C   G   T   N
         4, -2,  1, -2,  0,  // A (i=0)
        -2,  4, -2,  1,  0,  // C (i=1)
         1, -2,  4, -2,  0,  // G (i=2)
        -2,  1, -2,  4,  0,  // T (i=3)
         0,  0,  0,  0,  0   // N (i=4)
    };

    ProfileMatrix::ProfileMatrix() : len(0), dim(5), depth(0), prof() {}

    ProfileMatrix::ProfileMatrix(const std::string& seq)
        : len(static_cast<int>(seq.size())),
          dim(5),
          depth(seq.empty() ? 0 : 1),
          prof(static_cast<std::size_t>(len) * 5, 0U)
    {
        for (int i = 0; i < len; ++i) {
            const char ch = seq[static_cast<std::size_t>(i)];
            const int idx = baseIndex(ch);
            prof[static_cast<std::size_t>(i) * 5 + static_cast<std::size_t>(idx)] = 1U;
        }
    }

    ProfileMatrix::ProfileMatrix(const consensus::ConsensusJson& cj)
        : ProfileMatrix(fromConsensusCounts(cj))
    {
    }

    int ProfileMatrix::baseIndex(char ch)
    {
        switch (ch) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            case 'U': case 'u': return 3;
            case 'N': case 'n': return 4;
            default: return 4;
        }
    }

    bool ProfileMatrix::isGap(char ch)
    {
        return ch == '-' || ch == '.';
    }

    ProfileMatrix ProfileMatrix::fromAlignedSequence(const std::string& seq)
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

    ProfileMatrix ProfileMatrix::fromAlignedSequences(const std::vector<std::string>& aligned_sequences)
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

    ProfileMatrix ProfileMatrix::fromConsensusCounts(const consensus::ConsensusJson& cj)
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
            const std::uint64_t t_total = static_cast<std::uint64_t>(sc.t) +
                                          static_cast<std::uint64_t>(sc.u);
            pm.prof[off + 3] = static_cast<std::uint32_t>(t_total);
            pm.prof[off + 4] = sc.n;
        }

        return pm;
    }

    ProfileMatrix normalizeProfileEqualWeight(const ProfileMatrix& profile)
    {
        if (profile.len < 0 || profile.dim <= 0 ||
            profile.prof.size() != static_cast<std::size_t>(profile.len) * static_cast<std::size_t>(profile.dim)) {
            throw std::runtime_error("normalizeProfileEqualWeight: invalid profile shape");
        }

        ProfileMatrix normalized;
        normalized.len = profile.len;
        normalized.dim = profile.dim;
        normalized.depth = static_cast<int>(kProfileEqualWeightScale);
        normalized.prof.assign(profile.prof.size(), 0U);

        const std::uint64_t denom = profile.depth > 0
            ? static_cast<std::uint64_t>(profile.depth)
            : 1U;

        for (std::size_t i = 0; i < profile.prof.size(); ++i) {
            const std::uint64_t scaled =
                (static_cast<std::uint64_t>(profile.prof[i]) << kProfileEqualWeightShift) + (denom / 2U);
            const std::uint64_t value = scaled / denom;
            if (value > static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::runtime_error("normalizeProfileEqualWeight: normalized count overflow");
            }
            normalized.prof[i] = static_cast<std::uint32_t>(value);
        }

        return normalized;
    }

    ProfileMatrix combineProfilesEqualWeight(const std::vector<const ProfileMatrix*>& profiles)
    {
        if (profiles.empty()) {
            return ProfileMatrix();
        }

        const int len = profiles.front()->len;
        const int dim = profiles.front()->dim;
        if (len < 0 || dim <= 0) {
            throw std::runtime_error("combineProfilesEqualWeight: invalid profile shape");
        }

        ProfileMatrix combined;
        combined.len = len;
        combined.dim = dim;
        combined.depth = static_cast<int>(profiles.size() * kProfileEqualWeightScale);
        combined.prof.assign(static_cast<std::size_t>(len) * static_cast<std::size_t>(dim), 0U);

        for (const ProfileMatrix* profile : profiles) {
            if (profile == nullptr || profile->len != len || profile->dim != dim ||
                profile->prof.size() != static_cast<std::size_t>(len) * static_cast<std::size_t>(dim)) {
                throw std::runtime_error("combineProfilesEqualWeight: profile shape mismatch");
            }

            if (profile->depth == static_cast<int>(kProfileEqualWeightScale)) {
                for (std::size_t i = 0; i < profile->prof.size(); ++i) {
                    combined.prof[i] += profile->prof[i];
                }
                continue;
            }

            const std::uint64_t denom = profile->depth > 0
                ? static_cast<std::uint64_t>(profile->depth)
                : 1U;
            for (std::size_t i = 0; i < profile->prof.size(); ++i) {
                const std::uint64_t scaled =
                    (static_cast<std::uint64_t>(profile->prof[i]) << kProfileEqualWeightShift) + (denom / 2U);
                const std::uint64_t value = scaled / denom;
                if (value > static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max())) {
                    throw std::runtime_error("combineProfilesEqualWeight: normalized count overflow");
                }
                combined.prof[i] += static_cast<std::uint32_t>(value);
            }
        }

        return combined;
    }

    std::string profileToGappedSequence(const ProfileMatrix& profile)
    {
        if (profile.len < 0 || profile.dim < 5) {
            throw std::runtime_error("profileToGappedSequence: invalid profile shape");
        }

        const std::size_t len = static_cast<std::size_t>(profile.len);
        const std::size_t dim = static_cast<std::size_t>(profile.dim);
        if (profile.prof.size() < len * dim) {
            throw std::runtime_error("profileToGappedSequence: profile data is truncated");
        }

        static constexpr char bases[5] = {'A', 'C', 'G', 'T', 'N'};
        std::string seq;
        seq.reserve(len);

        for (std::size_t col = 0; col < len; ++col) {
            const std::size_t off = col * dim;
            std::uint64_t total = 0;
            std::uint32_t best_count = 0;
            std::size_t best_idx = 4;

            for (std::size_t idx = 0; idx < 5; ++idx) {
                const std::uint32_t count = profile.prof[off + idx];
                total += count;
                if (count > best_count) {
                    best_count = count;
                    best_idx = idx;
                }
            }

            seq.push_back(profile.depth - total > best_count ? '-' : bases[best_idx]);
        }

        return seq;
    }

    seq_io::SeqRecord reverseComplementRecord(const seq_io::SeqRecord& rec)
    {
        seq_io::SeqRecord out = rec;
        out.seq.resize(rec.seq.size());
        for (std::size_t i = 0; i < rec.seq.size(); ++i) {
            out.seq[i] = complementBase(rec.seq[rec.seq.size() - 1U - i]);
        }
        if (!out.qual.empty()) {
            std::reverse(out.qual.begin(), out.qual.end());
        }
        return out;
    }

    int auto_band(int qlen, int tlen, double indel_rate, int margin)
    {
        // 长度差异过大时不适合 banded DP
        if ((double)std::abs(qlen - tlen) / (double)std::max(qlen, tlen) > 0.5)
        {
            return -1;
        }
        // 经验公式：预期 indel 规模 + 安全边距
        return margin + static_cast<int>(indel_rate * (qlen + tlen / 2));
    }

    static int resolve_band_width(int configured_band_width, int ref_len, int query_len)
    {
        if (configured_band_width == AUTO_BAND_WIDTH) {
            return align::auto_band(ref_len, query_len);
        }
        return configured_band_width;
    }

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

        cfg.band_width = resolve_band_width(
            cfg.band_width,
            static_cast<int>(ref.size()),
            static_cast<int>(query.size()));

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

        cfg.band_width = resolve_band_width(
            cfg.band_width,
            ref.len,
            static_cast<int>(query.size()));

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

    cigar::Cigar_t globalAlignMM2(const std::string& ref,
                                  const std::string& query,
                                  const anchor::Anchors& anchors,
                                  align::AlignConfig cfg)
    {
        return globalAlignSeq2Seq(ref, query, anchors, cfg);
    }

    // 基于锚点的分段全局比对（minimap2 风格）
    // - 用锚点拆分为多个片段，逐段全局比对，最后合并
    cigar::Cigar_t globalAlignSeq2Seq(const std::string& ref,
                                      const std::string& query,
                                      const anchor::Anchors& anchors,
                                      align::AlignConfig cfg)
    {
        align::AlignConfig first_cfg = cfg;

        const std::size_t ref_len = ref.size();
        const std::size_t qry_len = query.size();

        // 链化锚点：获取最佳链
        anchor::Anchors sorted_anchors = anchors;
        anchor::ChainParams chain_params = anchor::default_chain_params();
        anchor::Anchors chain_anchors = anchor::chainAnchors(sorted_anchors, chain_params);
        if (chain_anchors.empty()) {
            return globalAlignKSW2(ref, query, cfg);
        }

        // 按 query 坐标排序
        std::sort(chain_anchors.begin(), chain_anchors.end(),
                  [](const anchor::Anchor& a, const anchor::Anchor& b) {
                      if (a.pos_qry != b.pos_qry) return a.pos_qry < b.pos_qry;
                      return a.pos_ref < b.pos_ref;
                  });

        const std::vector<MatchBlock> match_blocks = buildExtendedMatchBlocks(
            chain_anchors, ref_len, qry_len, query,
            [&ref](std::size_t pos) {
                return canonicalDnaBase(ref[pos]);
            });
        if (match_blocks.empty()) {
            return globalAlignKSW2(ref, query, cfg);
        }

        cigar::Cigar_t result;
        result.reserve(match_blocks.size() * 2 + 2);

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

        for (const MatchBlock& block : match_blocks) {
            append_segment(ref_pos, block.ref_start, qry_pos, block.qry_start, first_cfg);
            appendCigarOp(result, 'M', block.len);
            ref_pos = block.ref_start + block.len;
            qry_pos = block.qry_start + block.len;
        }

        // 右端：最后一个锚点到末尾
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg);

        // 最终一致性检查
        const std::size_t total_ref = cigar::getRefLength(result);
        const std::size_t total_qry = cigar::getQueryLength(result);
        if (total_ref != ref_len || total_qry != qry_len) {
            spdlog::error("globalAlignSeq2Seq: final cigar mismatch (ref:{}/{}, qry:{}/{}), fallback to global",
                         total_ref, ref_len, total_qry, qry_len);
            return globalAlignKSW2(ref, query, cfg);
        }

        return result;
    }

    cigar::Cigar_t globalAlignSeq2Profile(const ProfileMatrix& ref,
                                const std::string& ref_string,
                              const std::string& query,
                              const anchor::Anchors& anchors,
                              align::AlignConfig cfg)
    {
        align::AlignConfig first_cfg = cfg;

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

        const std::vector<MatchBlock> match_blocks = buildExtendedMatchBlocks(
            chain_anchors, ref_len, qry_len, query,
            [&ref](std::size_t pos) {
                return profileConsensusBaseAt(ref, pos);
            });
        if (match_blocks.empty()) {
            return globalAlignPSW(ref, query, cfg);
        }

        cigar::Cigar_t result;
        result.reserve(match_blocks.size() * 2 + 2);

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto append_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                  std::size_t qry_start, std::size_t qry_end,
                                  align::AlignConfig seg_cfg,
                                  bool use_terminal_ksw = false,
                                  bool reverse_from_anchor = false) {
            // 边界裁剪
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;

            const std::size_t seg_ref_len = ref_end - ref_start;
            const std::string seg_qry = query.substr(qry_start, qry_end - qry_start);
            const std::size_t seg_qry_len = seg_qry.size();

            cigar::Cigar_t seg_cigar;
            if (use_terminal_ksw) {
                seg_cigar = alignTerminalKSW2(
                    ref_string.substr(ref_start, seg_ref_len),
                    seg_qry,
                    seg_cfg,
                    reverse_from_anchor);
            } else {
                ProfileMatrix seg_ref = sliceProfile(ref, ref_start, ref_end);
                seg_cigar = globalAlignPSW(seg_ref, seg_qry, seg_cfg);
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
                    appendCigarOp(forced_cigar, 'I', seg_qry_len);
                }
                if (seg_ref_len > 0) {
                    appendCigarOp(forced_cigar, 'D', seg_ref_len);
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

        for (std::size_t i = 0; i < match_blocks.size(); ++i) {
            const MatchBlock& block = match_blocks[i];
            const bool left_terminal = (i == 0);
            append_segment(ref_pos, block.ref_start, qry_pos, block.qry_start,
                           left_terminal ? first_cfg : cfg,
                           left_terminal,
                           left_terminal);
            appendCigarOp(result, 'M', block.len);
            ref_pos = block.ref_start + block.len;
            qry_pos = block.qry_start + block.len;
        }

        // 右端：最后一个锚点到末尾
        append_segment(ref_pos, ref_len, qry_pos, qry_len, cfg, true, false);

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
                                               int thread,
                                               align::AlignConfig cfg)
    {
        align::AlignConfig first_cfg = cfg;

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

        const std::vector<MatchBlock> match_blocks = buildExtendedMatchBlocks(
            chain_anchors, ref_len, qry_len, query,
            [&ref](std::size_t pos) {
                return profileConsensusBaseAt(ref, pos);
            });
        if (match_blocks.empty()) {
            return globalAlignPSW(ref, query, cfg);
        }

        enum class SegmentKind {
            Psw,
            TerminalKsw,
            DirectMatch
        };

        struct SegmentTask {
            std::size_t ref_start = 0;
            std::size_t ref_end = 0;
            std::size_t qry_start = 0;
            std::size_t qry_end = 0;
            align::AlignConfig seg_cfg{};
            SegmentKind kind = SegmentKind::Psw;
            bool reverse_from_anchor = false;
        };

        std::vector<SegmentTask> tasks;
        tasks.reserve(match_blocks.size() * 2 + 2);

        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        auto push_segment = [&](std::size_t ref_start, std::size_t ref_end,
                                std::size_t qry_start, std::size_t qry_end,
                                align::AlignConfig seg_cfg,
                                SegmentKind kind,
                                bool reverse_from_anchor = false) {
            // 边界裁剪：保持和串行版 append_segment 完全一致
            ref_start = std::min(ref_start, ref_len);
            ref_end = std::min(ref_end, ref_len);
            qry_start = std::min(qry_start, qry_len);
            qry_end = std::min(qry_end, qry_len);

            if (ref_end < ref_start) ref_end = ref_start;
            if (qry_end < qry_start) qry_end = qry_start;
            if (ref_start == ref_end && qry_start == qry_end) {
                ref_pos = ref_end;
                qry_pos = qry_end;
                return;
            }

            SegmentTask task;
            task.ref_start = ref_start;
            task.ref_end = ref_end;
            task.qry_start = qry_start;
            task.qry_end = qry_end;
            task.seg_cfg = seg_cfg;
            task.kind = kind;
            task.reverse_from_anchor = reverse_from_anchor;

            tasks.push_back(std::move(task));

            // 串行版 append_segment 最终一定把游标推进到 ref_end/qry_end：
            // 1) 正常情况下 c_ref == seg_ref_len, c_qry == seg_qry_len
            // 2) fallback 情况下显式 ref_pos = ref_end, qry_pos = qry_end
            ref_pos = ref_end;
            qry_pos = qry_end;
        };

        auto push_match = [&](const MatchBlock& block) {
            SegmentTask task;
            task.ref_start = block.ref_start;
            task.ref_end = block.ref_start + block.len;
            task.qry_start = block.qry_start;
            task.qry_end = block.qry_start + block.len;
            task.seg_cfg = cfg;
            task.kind = SegmentKind::DirectMatch;
            tasks.push_back(std::move(task));
            ref_pos = block.ref_start + block.len;
            qry_pos = block.qry_start + block.len;
        };

        for (std::size_t i = 0; i < match_blocks.size(); ++i) {
            const MatchBlock& block = match_blocks[i];
            const bool left_terminal = (i == 0);
            push_segment(ref_pos, block.ref_start, qry_pos, block.qry_start,
                         left_terminal ? first_cfg : cfg,
                         left_terminal ? SegmentKind::TerminalKsw : SegmentKind::Psw,
                         left_terminal);
            push_match(block);
        }

        push_segment(ref_pos, ref_len, qry_pos, qry_len, cfg, SegmentKind::TerminalKsw, false);

        std::vector<cigar::Cigar_t> task_cigars(tasks.size());

    #ifdef _OPENMP
        const int use_threads = thread > 0 ? thread : omp_get_max_threads();

    #pragma omp parallel for default(none) shared(tasks, task_cigars, ref, ref_string, query) num_threads(use_threads) schedule(static)
    #endif
        for (int task_idx = 0; task_idx < static_cast<int>(tasks.size()); ++task_idx) {
            const SegmentTask& task = tasks[static_cast<std::size_t>(task_idx)];

            const std::size_t seg_ref_len = task.ref_end - task.ref_start;
            const std::size_t seg_qry_len = task.qry_end - task.qry_start;

            if (task.kind == SegmentKind::DirectMatch) {
                cigar::Cigar_t match_cigar;
                appendCigarOp(match_cigar, 'M', seg_ref_len);
                task_cigars[static_cast<std::size_t>(task_idx)] = std::move(match_cigar);
                continue;
            }

            std::string seg_qry = query.substr(task.qry_start, task.qry_end - task.qry_start);
            cigar::Cigar_t seg_cigar;
            if (task.kind == SegmentKind::TerminalKsw) {
                seg_cigar = alignTerminalKSW2(
                    ref_string.substr(task.ref_start, seg_ref_len),
                    seg_qry,
                    task.seg_cfg,
                    task.reverse_from_anchor);
            } else {
                ProfileMatrix seg_ref = sliceProfile(ref, task.ref_start, task.ref_end);
                seg_cigar = globalAlignPSW(seg_ref, seg_qry, task.seg_cfg);
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
                    appendCigarOp(forced_cigar, 'I', seg_qry_len);
                }

                if (seg_ref_len > 0) {
                    appendCigarOp(forced_cigar, 'D', seg_ref_len);
                }

                task_cigars[static_cast<std::size_t>(task_idx)] = std::move(forced_cigar);
            } else {
                task_cigars[static_cast<std::size_t>(task_idx)] = std::move(seg_cigar);
            }
        }

        // 按原始分段顺序合并，保证结果不受并行执行顺序影响
        cigar::Cigar_t result;
        result.reserve(match_blocks.size() * 2 + 2);

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
