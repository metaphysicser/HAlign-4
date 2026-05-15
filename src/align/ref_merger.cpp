#include "align.h"
#include "config.hpp"
#include "preprocess.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <omp.h>

namespace align {
namespace {

    struct InsertionSegment {
        std::size_t slot = 0;
        std::string bases;
    };

    struct InsertionKey {
        std::size_t record_index = 0;
        std::size_t slot = 0;

        bool operator==(const InsertionKey& other) const noexcept
        {
            return record_index == other.record_index && slot == other.slot;
        }
    };

    struct InsertionKeyHash {
        std::size_t operator()(const InsertionKey& key) const noexcept
        {
            std::size_t seed = std::hash<std::size_t>{}(key.record_index);
            const std::size_t value = std::hash<std::size_t>{}(key.slot);
            seed ^= value + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) +
                    (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };

    struct SlotInsertionEntry {
        InsertionKey key;
        std::string id;
        std::string bases;
    };

    struct ReferenceGuidedInsertionMsa {
        std::vector<std::size_t> widths;
        std::unordered_map<InsertionKey, std::string, InsertionKeyHash> aligned_segments;
        std::size_t records = 0;
        std::size_t segments = 0;
        std::size_t aligned_slots = 0;
    };

    struct WriteState {
        std::size_t expected_length = 0;
        bool length_initialized = false;
        std::size_t seq_count = 0;
    };

    struct CoordinateTransform {
        std::size_t local_length = 0;
        std::size_t final_length = 0;
        bool identity = true;
        cigar::Cigar_t local_to_final_cigar;
        std::vector<std::size_t> slot_to_final;

        std::size_t mapSlot(std::size_t local_slot) const
        {
            if (identity) {
                if (local_slot > final_length) {
                    throw std::runtime_error("CoordinateTransform: local slot out of range");
                }
                return local_slot;
            }
            if (local_slot >= slot_to_final.size()) {
                throw std::runtime_error("CoordinateTransform: local slot out of range");
            }
            return slot_to_final[local_slot];
        }
    };

    struct ReferenceContext {
        std::string consensus_id;
        std::string base_reference_seq;
        bool profile_alignment_mode = false;
        CoordinateTransform identity_transform;
        CoordinateTransform consensus_transform;
        std::unordered_map<std::string, CoordinateTransform> ref_transforms;

        std::size_t baseLength() const noexcept
        {
            return base_reference_seq.size();
        }

        const CoordinateTransform& transformFor(const std::string& ref_name) const
        {
            if (ref_name == consensus_id) {
                return profile_alignment_mode ? identity_transform : consensus_transform;
            }

            const auto it = ref_transforms.find(ref_name);
            if (it != ref_transforms.end()) {
                return it->second;
            }

            throw std::runtime_error("Reference '" + ref_name + "' not found in reference coordinate map");
        }
    };

    void appendCigarRun(cigar::Cigar_t& out, char op, std::uint32_t len)
    {
        if (len == 0) {
            return;
        }
        cigar::appendCigar(out, cigar::Cigar_t{cigar::cigarToInt(op, len)});
    }

    cigar::Cigar_t cigarFromAlignedSequence(std::string_view aligned_seq)
    {
        cigar::Cigar_t result;
        if (aligned_seq.empty()) {
            return result;
        }

        char current_op = '\0';
        std::uint32_t current_len = 0;

        for (const char base : aligned_seq) {
            const char op = (base == '-') ? 'D' : 'M';
            if (op == current_op) {
                ++current_len;
                continue;
            }

            appendCigarRun(result, current_op, current_len);
            current_op = op;
            current_len = 1;
        }
        appendCigarRun(result, current_op, current_len);

        return result;
    }

    CoordinateTransform makeIdentityTransform(std::size_t len)
    {
        CoordinateTransform transform;
        transform.local_length = len;
        transform.final_length = len;
        transform.identity = true;
        return transform;
    }

    CoordinateTransform makeTransformFromCigar(cigar::Cigar_t cigar)
    {
        CoordinateTransform transform;
        transform.local_length = cigar::getQueryLength(cigar);
        transform.final_length = cigar::getRefLength(cigar);
        transform.identity = false;
        transform.local_to_final_cigar = std::move(cigar);
        transform.slot_to_final.assign(transform.local_length + 1U, 0U);

        std::size_t local_pos = 0;
        std::size_t final_pos = 0;
        transform.slot_to_final[0] = 0;

        for (const cigar::CigarUnit unit : transform.local_to_final_cigar) {
            char op = '\0';
            std::uint32_t len = 0;
            cigar::intToCigar(unit, op, len);

            switch (op) {
                case 'M':
                case '=':
                case 'X':
                    for (std::uint32_t i = 0; i < len; ++i) {
                        transform.slot_to_final[local_pos] = final_pos;
                        ++local_pos;
                        ++final_pos;
                        transform.slot_to_final[local_pos] = final_pos;
                    }
                    break;
                case 'D':
                case 'N':
                    final_pos += static_cast<std::size_t>(len);
                    transform.slot_to_final[local_pos] = final_pos;
                    break;
                case 'H':
                case 'P':
                    break;
                default:
                    throw std::runtime_error("Coordinate transform CIGAR contains unsupported op: " +
                                             std::string(1, op));
            }
        }

        if (local_pos != transform.local_length || final_pos != transform.final_length) {
            throw std::runtime_error("Coordinate transform length mismatch");
        }

        return transform;
    }

    CoordinateTransform makeTransformFromAlignedSequence(std::string_view aligned_seq)
    {
        return makeTransformFromCigar(cigarFromAlignedSequence(aligned_seq));
    }

    bool hasUsableCigar(const seq_io::SamRecord& sam_rec)
    {
        return !sam_rec.cigar.empty() && sam_rec.cigar != "*";
    }

    bool samPathHasRecords(const FilePath& sam_path)
    {
        if (!std::filesystem::exists(sam_path)) {
            throw std::runtime_error("SAM file does not exist: " + sam_path.string());
        }
        return std::filesystem::file_size(sam_path) != 0;
    }

    std::vector<FilePath> concatPaths(const std::vector<FilePath>& first,
                                      const std::vector<FilePath>& second)
    {
        std::vector<FilePath> out;
        out.reserve(first.size() + second.size());
        out.insert(out.end(), first.begin(), first.end());
        out.insert(out.end(), second.begin(), second.end());
        return out;
    }

    void applyTransformToBaseCoordinate(std::string& seq, const CoordinateTransform& transform)
    {
        if (transform.identity) {
            if (seq.size() != transform.final_length) {
                throw std::runtime_error("Projected sequence length " + std::to_string(seq.size()) +
                                         " does not match reference length " +
                                         std::to_string(transform.final_length));
            }
            return;
        }

        cigar::padQueryToRefByCigar(seq, transform.local_to_final_cigar);
        if (seq.size() != transform.final_length) {
            throw std::runtime_error("Coordinate transform produced length " + std::to_string(seq.size()) +
                                     ", expected " + std::to_string(transform.final_length));
        }
    }

    seq_io::SeqRecord projectSamToBaseReference(const seq_io::SamRecord& sam_rec,
                                                const ReferenceContext& context)
    {
        seq_io::SeqRecord fasta_rec = seq_io::samRecordToSeqRecord(sam_rec, false);
        const CoordinateTransform& transform = context.transformFor(sam_rec.rname);

        if (hasUsableCigar(sam_rec)) {
            cigar::Cigar_t cigar_ops = cigar::stringToCigar(sam_rec.cigar);
            cigar::delQueryToRefByCigar(fasta_rec.seq, cigar_ops);
        }

        applyTransformToBaseCoordinate(fasta_rec.seq, transform);
        return fasta_rec;
    }

    template <typename Segment>
    void sortAndMergeBySlot(std::vector<Segment>& segments)
    {
        std::sort(segments.begin(), segments.end(),
                  [](const Segment& a, const Segment& b) {
                      return a.slot < b.slot;
                  });

        std::size_t write_pos = 0;
        for (std::size_t read_pos = 0; read_pos < segments.size(); ++read_pos) {
            if (write_pos != 0 && segments[write_pos - 1].slot == segments[read_pos].slot) {
                if constexpr (std::is_same_v<Segment, InsertionSegment>) {
                    segments[write_pos - 1].bases += segments[read_pos].bases;
                } else {
                    segments[write_pos - 1].len += segments[read_pos].len;
                }
                continue;
            }
            if (write_pos != read_pos) {
                segments[write_pos] = std::move(segments[read_pos]);
            }
            ++write_pos;
        }
        segments.resize(write_pos);
    }

    std::vector<InsertionSegment> collectInsertionSegments(const seq_io::SamRecord& sam_rec,
                                                           const CoordinateTransform& transform)
    {
        std::vector<InsertionSegment> segments;
        if (!hasUsableCigar(sam_rec)) {
            return segments;
        }

        const cigar::Cigar_t cigar_ops = cigar::stringToCigar(sam_rec.cigar);
        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        for (const cigar::CigarUnit unit : cigar_ops) {
            char op = '\0';
            std::uint32_t len = 0;
            cigar::intToCigar(unit, op, len);

            switch (op) {
                case 'M':
                case '=':
                case 'X':
                    ref_pos += static_cast<std::size_t>(len);
                    qry_pos += static_cast<std::size_t>(len);
                    break;
                case 'D':
                case 'N':
                    ref_pos += static_cast<std::size_t>(len);
                    break;
                case 'I': {
                    const std::size_t insert_len = static_cast<std::size_t>(len);
                    if (qry_pos + insert_len > sam_rec.seq.size()) {
                        throw std::runtime_error("CIGAR insertion consumes past query end for " + sam_rec.qname);
                    }
                    segments.push_back(InsertionSegment{
                        transform.mapSlot(ref_pos),
                        sam_rec.seq.substr(qry_pos, insert_len)
                    });
                    qry_pos += insert_len;
                    break;
                }
                case 'S':
                    qry_pos += static_cast<std::size_t>(len);
                    break;
                case 'H':
                case 'P':
                    break;
                default:
                    throw std::runtime_error("Unsupported CIGAR op: " + std::string(1, op));
            }
        }

        sortAndMergeBySlot(segments);
        return segments;
    }

    void appendInsertionSlot(std::string& out,
                             const std::vector<InsertionSegment>& segments,
                             std::size_t& segment_idx,
                             std::size_t slot,
                             std::size_t width)
    {
        std::size_t observed = 0;
        if (segment_idx < segments.size() && segments[segment_idx].slot == slot) {
            observed = segments[segment_idx].bases.size();
            if (observed > width) {
                throw std::runtime_error("Insertion segment length exceeds allocated slot width");
            }
            out.append(segments[segment_idx].bases);
            ++segment_idx;
        }
        out.append(width - observed, '-');
    }

    std::size_t expandedLength(std::size_t base_length, const std::vector<std::size_t>& insertion_widths)
    {
        std::size_t len = base_length;
        for (const std::size_t width : insertion_widths) {
            len += width;
        }
        return len;
    }

    std::string expandBaseWithInsertions(const std::string& base_seq,
                                         const std::vector<InsertionSegment>& segments,
                                         const std::vector<std::size_t>& insertion_widths)
    {
        if (insertion_widths.size() != base_seq.size() + 1U) {
            throw std::runtime_error("Insertion width vector does not match base reference length");
        }

        std::string out;
        out.reserve(expandedLength(base_seq.size(), insertion_widths));

        std::size_t segment_idx = 0;
        appendInsertionSlot(out, segments, segment_idx, 0, insertion_widths[0]);
        for (std::size_t i = 0; i < base_seq.size(); ++i) {
            out.push_back(base_seq[i]);
            appendInsertionSlot(out, segments, segment_idx, i + 1U, insertion_widths[i + 1U]);
        }

        if (segment_idx != segments.size()) {
            throw std::runtime_error("Insertion segment mapped outside the reference coordinate");
        }

        return out;
    }

    std::string expandBaseWithInsertionGaps(const std::string& base_seq,
                                            const std::vector<std::size_t>& insertion_widths)
    {
        static const std::vector<InsertionSegment> empty_segments;
        return expandBaseWithInsertions(base_seq, empty_segments, insertion_widths);
    }

    void writeChecked(seq_io::SeqWriter& writer,
                      const seq_io::SeqRecord& rec,
                      WriteState& state,
                      ProgressBar& progress)
    {
        if (!state.length_initialized) {
            state.expected_length = rec.seq.size();
            state.length_initialized = true;
        } else if (rec.seq.size() != state.expected_length) {
            throw std::runtime_error("Sequence length mismatch: " + rec.id +
                                     " length " + std::to_string(rec.seq.size()) +
                                     ", expected " + std::to_string(state.expected_length));
        }

        writer.writeFasta(rec);
        ++state.seq_count;
        progress.tick();
    }

    std::size_t writeReferenceRecords(seq_io::SeqWriter& writer,
                                      const FilePath& aligned_reference_fasta,
                                      const ReferenceContext& context,
                                      const std::vector<std::size_t>& insertion_widths,
                                      WriteState& state,
                                      ProgressBar& progress)
    {
        seq_io::KseqReader reader(aligned_reference_fasta);
        seq_io::SeqRecord rec;
        std::size_t written = 0;

        while (reader.next(rec)) {
            seq_io::cleanSequence(rec);
            if (rec.seq.size() != context.baseLength()) {
                throw std::runtime_error("Reference coordinate length mismatch for " + rec.id +
                                         ": " + std::to_string(rec.seq.size()) +
                                         " vs " + std::to_string(context.baseLength()));
            }
            rec.seq = expandBaseWithInsertionGaps(rec.seq, insertion_widths);
            writeChecked(writer, rec, state, progress);
            ++written;
        }

        if (written == 0) {
            throw std::runtime_error("No reference records found in " + aligned_reference_fasta.string());
        }

        writer.flush();
        return written;
    }

    std::vector<std::size_t> zeroInsertionWidths(std::size_t base_length)
    {
        return std::vector<std::size_t>(base_length + 1U, 0U);
    }

    std::string makeInsertionMsaId(const InsertionKey& key)
    {
        return "rgins_" + std::to_string(key.record_index) + "_" +
               std::to_string(key.slot);
    }

    std::size_t parseSizeToken(std::string_view token, std::string_view full_id)
    {
        if (token.empty()) {
            throw std::runtime_error("Invalid reference-guided insertion MSA id: " +
                                     std::string(full_id));
        }

        std::size_t consumed = 0;
        unsigned long long value = 0;
        try {
            value = std::stoull(std::string(token), &consumed);
        } catch (...) {
            throw std::runtime_error("Invalid reference-guided insertion MSA id: " +
                                     std::string(full_id));
        }

        if (consumed != token.size() ||
            value > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
            throw std::runtime_error("Invalid reference-guided insertion MSA id: " +
                                     std::string(full_id));
        }

        return static_cast<std::size_t>(value);
    }

    InsertionKey parseInsertionMsaId(std::string_view id)
    {
        constexpr std::string_view prefix = "rgins_";
        if (id.substr(0, prefix.size()) != prefix) {
            throw std::runtime_error("Unexpected reference-guided insertion MSA id: " +
                                     std::string(id));
        }

        const std::size_t slot_sep = id.find('_', prefix.size());
        if (slot_sep == std::string_view::npos) {
            throw std::runtime_error("Invalid reference-guided insertion MSA id: " +
                                     std::string(id));
        }

        return InsertionKey{
            parseSizeToken(id.substr(prefix.size(), slot_sep - prefix.size()), id),
            parseSizeToken(id.substr(slot_sep + 1U), id)
        };
    }

    std::string normalizedAlignedInsertion(std::string seq)
    {
        for (char& ch : seq) {
            if (ch == '.') {
                ch = '-';
            } else {
                ch = static_cast<char>(seq_io::clean_table[static_cast<unsigned char>(ch)]);
            }
        }
        return seq;
    }

    std::string ungappedInsertion(std::string_view seq)
    {
        std::string out;
        out.reserve(seq.size());
        for (const char ch : seq) {
            if (ch == '-' || ch == '.') {
                continue;
            }
            out.push_back(static_cast<char>(
                seq_io::clean_table[static_cast<unsigned char>(ch)]));
        }
        return out;
    }

    void writeSlotInsertionMsaInput(const FilePath& input_path,
                                    const std::vector<SlotInsertionEntry>& entries)
    {
        seq_io::SeqWriter writer(input_path, 80);
        for (const SlotInsertionEntry& entry : entries) {
            seq_io::SeqRecord rec;
            rec.id = entry.id;
            rec.seq = entry.bases;
            writer.writeFasta(rec);
        }
        writer.flush();
    }

    std::size_t readSlotInsertionMsaOutput(
        const FilePath& aligned_path,
        const std::vector<SlotInsertionEntry>& entries,
        std::unordered_map<InsertionKey, std::string, InsertionKeyHash>& aligned_segments)
    {
        std::unordered_map<std::string, std::size_t> expected;
        expected.reserve(entries.size());
        for (std::size_t i = 0; i < entries.size(); ++i) {
            expected.emplace(entries[i].id, i);
        }

        std::vector<bool> seen(entries.size(), false);
        seq_io::KseqReader reader(aligned_path);
        seq_io::SeqRecord rec;
        std::size_t read_count = 0;
        std::size_t width = 0;
        bool width_initialized = false;

        while (reader.next(rec)) {
            const InsertionKey parsed_key = parseInsertionMsaId(rec.id);
            const auto it = expected.find(rec.id);
            if (it == expected.end()) {
                throw std::runtime_error(
                    "External insertion MSA returned an unknown sequence id: " + rec.id);
            }

            const std::size_t entry_idx = it->second;
            const SlotInsertionEntry& entry = entries[entry_idx];
            if (!(parsed_key == entry.key)) {
                throw std::runtime_error(
                    "External insertion MSA sequence id does not match its slot metadata: " +
                    rec.id);
            }
            if (seen[entry_idx]) {
                throw std::runtime_error(
                    "External insertion MSA returned duplicate sequence id: " + rec.id);
            }
            seen[entry_idx] = true;

            std::string aligned_seq = normalizedAlignedInsertion(std::move(rec.seq));
            if (ungappedInsertion(aligned_seq) != ungappedInsertion(entry.bases)) {
                throw std::runtime_error(
                    "External insertion MSA changed insertion sequence: " + rec.id);
            }

            if (!width_initialized) {
                width = aligned_seq.size();
                width_initialized = true;
            } else if (aligned_seq.size() != width) {
                throw std::runtime_error(
                    "External insertion MSA produced inconsistent lengths for slot " +
                    std::to_string(entry.key.slot));
            }

            aligned_segments.emplace(entry.key, std::move(aligned_seq));
            ++read_count;
        }

        if (read_count != entries.size()) {
            throw std::runtime_error(
                "External insertion MSA returned " + std::to_string(read_count) +
                " records, expected " + std::to_string(entries.size()));
        }

        return width;
    }

    ReferenceGuidedInsertionMsa alignReferenceGuidedInsertionSlots(
        const std::vector<FilePath>& insertion_sam_paths,
        const ReferenceContext& context,
        const FilePath& msa_dir,
        const std::string& msa_tool,
        int threads)
    {
        ReferenceGuidedInsertionMsa result;
        result.widths = zeroInsertionWidths(context.baseLength());

        std::unordered_map<std::size_t, std::vector<SlotInsertionEntry>> slot_groups;

        for (const FilePath& sam_path : insertion_sam_paths) {
            if (!samPathHasRecords(sam_path)) {
                continue;
            }

            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;
            while (reader.next(sam_rec)) {
                const std::size_t record_index = result.records++;
                const CoordinateTransform& transform = context.transformFor(sam_rec.rname);
                const std::vector<InsertionSegment> segments =
                    collectInsertionSegments(sam_rec, transform);

                for (const InsertionSegment& segment : segments) {
                    if (segment.slot >= result.widths.size()) {
                        throw std::runtime_error("Insertion slot outside reference coordinate");
                    }

                    InsertionKey key{record_index, segment.slot};
                    SlotInsertionEntry entry;
                    entry.key = key;
                    entry.id = makeInsertionMsaId(key);
                    entry.bases = segment.bases;
                    slot_groups[segment.slot].push_back(std::move(entry));
                    ++result.segments;
                }
            }
        }

        if (result.segments == 0) {
            return result;
        }

        result.aligned_segments.reserve(result.segments);

        if (std::filesystem::exists(msa_dir)) {
            file_io::removeAll(msa_dir);
        }
        file_io::ensureDirectoryExists(msa_dir, "reference-guided insertion MSA dir");

        const std::string effective_msa_tool = msa_tool.empty() ? DEFAULT_MSA_CMD : msa_tool;

        for (const auto& [slot, entries] : slot_groups) {
            if (entries.empty()) {
                continue;
            }

            if (entries.size() == 1U) {
                const SlotInsertionEntry& entry = entries.front();
                const std::string aligned_seq = normalizedAlignedInsertion(entry.bases);
                result.widths[slot] = aligned_seq.size();
                result.aligned_segments.emplace(entry.key, aligned_seq);
                continue;
            }

            const FilePath input_path =
                msa_dir / ("slot_" + std::to_string(slot) + ".fasta");
            const FilePath aligned_path =
                msa_dir / ("slot_" + std::to_string(slot) + "_aligned.fasta");

            writeSlotInsertionMsaInput(input_path, entries);
            alignConsensusSequence(input_path, aligned_path, effective_msa_tool, threads);
            file_io::requireRegularFile(aligned_path, "aligned insertion MSA");

            result.widths[slot] =
                readSlotInsertionMsaOutput(aligned_path, entries, result.aligned_segments);
            ++result.aligned_slots;
        }

        return result;
    }

    std::size_t writeInsertionTsv(const std::vector<FilePath>& insertion_sam_paths,
                                  const ReferenceContext& context,
                                  const FilePath& insertion_tsv_path)
    {
        if (insertion_tsv_path.empty()) {
            return 0;
        }

        file_io::ensureParentDirExists(insertion_tsv_path);
        std::ofstream out(insertion_tsv_path, std::ios::out | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("failed to open insertion TSV: " + insertion_tsv_path.string());
        }

        out << "query_id\tref_id\tref_slot\tref_pos_1based_after\tinsertion_length\tinsertion_sequence\n";

        std::size_t rows = 0;
        for (const FilePath& sam_path : insertion_sam_paths) {
            if (!samPathHasRecords(sam_path)) {
                continue;
            }

            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;
            while (reader.next(sam_rec)) {
                const CoordinateTransform& transform = context.transformFor(sam_rec.rname);
                const std::vector<InsertionSegment> segments =
                    collectInsertionSegments(sam_rec, transform);
                for (const InsertionSegment& segment : segments) {
                    out << sam_rec.qname << '\t'
                        << sam_rec.rname << '\t'
                        << segment.slot << '\t'
                        << segment.slot << '\t'
                        << segment.bases.size() << '\t'
                        << segment.bases << '\n';
                    ++rows;
                }
            }
        }

        out.flush();
        if (!out) {
            throw std::runtime_error("failed to write insertion TSV: " + insertion_tsv_path.string());
        }

        spdlog::info("Insertion TSV written: {} rows, {}", rows, insertion_tsv_path.string());
        return rows;
    }

    seq_io::SeqRecord projectSamRecord(const seq_io::SamRecord& sam_rec,
                                       const ReferenceContext& context,
                                       const std::vector<std::size_t>& insertion_widths,
                                       bool include_record_insertions)
    {
        seq_io::SeqRecord rec = projectSamToBaseReference(sam_rec, context);
        const CoordinateTransform& transform = context.transformFor(sam_rec.rname);

        if (include_record_insertions) {
            rec.seq = expandBaseWithInsertions(
                rec.seq, collectInsertionSegments(sam_rec, transform), insertion_widths);
        } else {
            rec.seq = expandBaseWithInsertionGaps(rec.seq, insertion_widths);
        }

        return rec;
    }

    std::vector<InsertionSegment> lookupAlignedInsertionSegments(
        const seq_io::SamRecord& sam_rec,
        const CoordinateTransform& transform,
        std::size_t record_index,
        const std::vector<std::size_t>& insertion_widths,
        const std::unordered_map<InsertionKey, std::string, InsertionKeyHash>& aligned_segments)
    {
        std::vector<InsertionSegment> raw_segments =
            collectInsertionSegments(sam_rec, transform);

        for (InsertionSegment& segment : raw_segments) {
            InsertionKey key{record_index, segment.slot};
            const auto it = aligned_segments.find(key);
            if (it == aligned_segments.end()) {
                throw std::runtime_error(
                    "Aligned insertion segment missing for record " +
                    std::to_string(record_index) + ", slot " +
                    std::to_string(segment.slot));
            }
            if (segment.slot >= insertion_widths.size()) {
                throw std::runtime_error("Insertion slot outside reference coordinate");
            }
            if (it->second.size() != insertion_widths[segment.slot]) {
                throw std::runtime_error(
                    "Aligned insertion segment width mismatch for record " +
                    std::to_string(record_index) + ", slot " +
                    std::to_string(segment.slot));
            }
            segment.bases = it->second;
        }

        return raw_segments;
    }

    seq_io::SeqRecord projectSamRecordWithAlignedInsertions(
        const seq_io::SamRecord& sam_rec,
        const ReferenceContext& context,
        const std::vector<std::size_t>& insertion_widths,
        std::size_t record_index,
        const std::unordered_map<InsertionKey, std::string, InsertionKeyHash>& aligned_segments)
    {
        seq_io::SeqRecord rec = projectSamToBaseReference(sam_rec, context);
        const CoordinateTransform& transform = context.transformFor(sam_rec.rname);
        rec.seq = expandBaseWithInsertions(
            rec.seq,
            lookupAlignedInsertionSegments(
                sam_rec, transform, record_index, insertion_widths, aligned_segments),
            insertion_widths);
        return rec;
    }

    void writeProjectedSamPaths(
        const std::vector<FilePath>& sam_paths,
        const ReferenceContext& context,
        const std::vector<std::size_t>& insertion_widths,
        bool include_record_insertions,
        std::size_t batch_size,
        int threads,
        seq_io::SeqWriter& writer,
        WriteState& state,
        ProgressBar& progress)
    {
        constexpr std::size_t default_batch_size = 2560;
        const std::size_t effective_batch_size = batch_size > 0 ? batch_size : default_batch_size;
        const int use_threads = std::max(1, threads);

        std::vector<seq_io::SamRecord> sam_batch;
        std::vector<seq_io::SeqRecord> fasta_batch;
        std::vector<std::string> errors;
        sam_batch.reserve(effective_batch_size);

        for (const FilePath& sam_path : sam_paths) {
            if (!samPathHasRecords(sam_path)) {
                continue;
            }

            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;

            while (true) {
                sam_batch.clear();
                while (sam_batch.size() < effective_batch_size && reader.next(sam_rec)) {
                    sam_batch.push_back(sam_rec);
                }
                if (sam_batch.empty()) {
                    break;
                }

                const std::size_t current_size = sam_batch.size();
                fasta_batch.clear();
                fasta_batch.resize(current_size);
                errors.clear();
                errors.resize(current_size);

#pragma omp parallel for default(none) shared(sam_batch, fasta_batch, errors, context, insertion_widths, current_size, include_record_insertions) num_threads(use_threads) schedule(dynamic, 8)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(current_size); ++i) {
                    try {
                        fasta_batch[static_cast<std::size_t>(i)] = projectSamRecord(
                            sam_batch[static_cast<std::size_t>(i)],
                            context,
                            insertion_widths,
                            include_record_insertions);
                    } catch (const std::exception& e) {
                        errors[static_cast<std::size_t>(i)] = e.what();
                    } catch (...) {
                        errors[static_cast<std::size_t>(i)] = "unknown projection error";
                    }
                }

                for (const std::string& err : errors) {
                    if (!err.empty()) {
                        throw std::runtime_error(err);
                    }
                }

                for (const seq_io::SeqRecord& rec : fasta_batch) {
                    writeChecked(writer, rec, state, progress);
                }
                writer.flush();
            }
        }
    }

    void writeProjectedSamPathsWithAlignedInsertions(
        const std::vector<FilePath>& sam_paths,
        const ReferenceContext& context,
        const ReferenceGuidedInsertionMsa& insertion_msa,
        std::size_t batch_size,
        int threads,
        seq_io::SeqWriter& writer,
        WriteState& state,
        ProgressBar& progress)
    {
        constexpr std::size_t default_batch_size = 2560;
        const std::size_t effective_batch_size = batch_size > 0 ? batch_size : default_batch_size;
        const int use_threads = std::max(1, threads);

        std::vector<seq_io::SamRecord> sam_batch;
        std::vector<std::size_t> record_indices;
        std::vector<seq_io::SeqRecord> fasta_batch;
        std::vector<std::string> errors;
        sam_batch.reserve(effective_batch_size);
        record_indices.reserve(effective_batch_size);

        std::size_t next_record_index = 0;

        for (const FilePath& sam_path : sam_paths) {
            if (!samPathHasRecords(sam_path)) {
                continue;
            }

            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;

            while (true) {
                sam_batch.clear();
                record_indices.clear();
                while (sam_batch.size() < effective_batch_size && reader.next(sam_rec)) {
                    record_indices.push_back(next_record_index++);
                    sam_batch.push_back(sam_rec);
                }
                if (sam_batch.empty()) {
                    break;
                }

                const std::size_t current_size = sam_batch.size();
                fasta_batch.clear();
                fasta_batch.resize(current_size);
                errors.clear();
                errors.resize(current_size);

#pragma omp parallel for default(none) shared(sam_batch, record_indices, fasta_batch, errors, context, insertion_msa, current_size) num_threads(use_threads) schedule(dynamic, 8)
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(current_size); ++i) {
                    try {
                        fasta_batch[static_cast<std::size_t>(i)] =
                            projectSamRecordWithAlignedInsertions(
                                sam_batch[static_cast<std::size_t>(i)],
                                context,
                                insertion_msa.widths,
                                record_indices[static_cast<std::size_t>(i)],
                                insertion_msa.aligned_segments);
                    } catch (const std::exception& e) {
                        errors[static_cast<std::size_t>(i)] = e.what();
                    } catch (...) {
                        errors[static_cast<std::size_t>(i)] = "unknown projection error";
                    }
                }

                for (const std::string& err : errors) {
                    if (!err.empty()) {
                        throw std::runtime_error(err);
                    }
                }

                for (const seq_io::SeqRecord& rec : fasta_batch) {
                    writeChecked(writer, rec, state, progress);
                }
                writer.flush();
            }
        }

        if (next_record_index != insertion_msa.records) {
            throw std::runtime_error(
                "Reference-guided insertion record count changed between MSA and output");
        }
    }

    std::size_t writeExternalMsaInput(const std::vector<FilePath>& insertion_sam_paths,
                                      const ReferenceContext& context,
                                      const FilePath& insertion_fasta_path)
    {
        seq_io::SeqWriter writer(insertion_fasta_path, 80);

        seq_io::SeqRecord reference_rec;
        reference_rec.id = context.consensus_id;
        reference_rec.seq = context.base_reference_seq;
        writer.writeFasta(reference_rec);

        std::size_t query_count = 0;
        for (const FilePath& sam_path : insertion_sam_paths) {
            if (!samPathHasRecords(sam_path)) {
                continue;
            }

            seq_io::SamReader reader(sam_path);
            seq_io::SamRecord sam_rec;
            while (reader.next(sam_rec)) {
                writer.writeFasta(seq_io::samRecordToSeqRecord(sam_rec, false));
                ++query_count;
            }
        }

        writer.flush();
        return query_count;
    }

    seq_io::SeqRecord readFirstFastaRecord(const FilePath& fasta_path)
    {
        seq_io::KseqReader reader(fasta_path);
        seq_io::SeqRecord rec;
        if (!reader.next(rec)) {
            throw std::runtime_error("FASTA is empty: " + fasta_path.string());
        }
        seq_io::cleanSequence(rec);
        return rec;
    }

    std::vector<std::size_t> insertionWidthsFromAlignedReference(
        const std::string& base_reference_seq,
        const std::string& aligned_reference_seq)
    {
        std::vector<std::size_t> widths(base_reference_seq.size() + 1U, 0U);
        std::size_t base_pos = 0;

        for (const char aligned_base : aligned_reference_seq) {
            if (base_pos < base_reference_seq.size() &&
                aligned_base == base_reference_seq[base_pos]) {
                ++base_pos;
                continue;
            }

            if (aligned_base == '-') {
                ++widths[base_pos];
                continue;
            }

            throw std::runtime_error("External insertion MSA changed the reference sequence");
        }

        if (base_pos != base_reference_seq.size()) {
            throw std::runtime_error("External insertion MSA reference is truncated");
        }

        return widths;
    }

    std::size_t writeAlignedInsertionRecords(seq_io::SeqWriter& writer,
                                             const FilePath& aligned_insertion_fasta,
                                             WriteState& state,
                                             ProgressBar& progress)
    {
        seq_io::KseqReader reader(aligned_insertion_fasta);
        seq_io::SeqRecord rec;
        bool skip_reference = true;
        std::size_t written = 0;

        while (reader.next(rec)) {
            if (skip_reference) {
                skip_reference = false;
                continue;
            }
            seq_io::cleanSequence(rec);
            writeChecked(writer, rec, state, progress);
            ++written;
        }

        writer.flush();
        return written;
    }

} // namespace

    void RefAligner::parseAlignedReferencesToCigar(
        const FilePath& aligned_fasta_path,
        std::unordered_map<std::string, cigar::Cigar_t>& out_ref_aligned_map,
        std::vector<bool>& out_ref_gap_pos) const
    {
        out_ref_aligned_map.clear();
        out_ref_gap_pos.clear();

        seq_io::KseqReader reader(aligned_fasta_path);
        seq_io::SeqRecord rec;
        std::size_t parsed_seq_count = 0;

        while (reader.next(rec)) {
            if (parsed_seq_count == 0) {
                out_ref_gap_pos.reserve(rec.seq.size());
                for (const char base : rec.seq) {
                    out_ref_gap_pos.push_back(base == '-');
                }
            }

            out_ref_aligned_map[rec.id] = cigarFromAlignedSequence(rec.seq);
            ++parsed_seq_count;
        }

        if (parsed_seq_count == 0) {
            throw std::runtime_error(
                "parseAlignedReferencesToCigar: input FASTA is empty: " +
                aligned_fasta_path.string());
        }
    }

    void RefAligner::mergeAlignedResults(const FilePath output, const MergeOptions& options)
    {
        ProgressBar progress("merge");

        const FilePath result_dir = work_dir / RESULTS_DIR;
        const FilePath consensus_aligned_file =
            FilePath(work_dir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath insertion_fasta_path = result_dir / ALL_INSERTION_FASTA;
        const FilePath aligned_insertion_fasta = result_dir / ALIGNED_INSERTION_FASTA;
        const FilePath reference_guided_insertion_dir =
            result_dir / "reference_guided_insertions";

        std::unordered_map<std::string, cigar::Cigar_t> ref_aligned_map;
        std::vector<bool> ref_gap_pos;
        parseAlignedReferencesToCigar(consensus_aligned_file, ref_aligned_map, ref_gap_pos);

        ReferenceContext context;
        context.consensus_id = consensus_seq.id;
        context.profile_alignment_mode = profile_alignment_mode;
        context.base_reference_seq = consensus_gap_seq.seq.empty() ? consensus_seq.seq : consensus_gap_seq.seq;
        context.identity_transform = makeIdentityTransform(context.baseLength());
        context.consensus_transform = profile_alignment_mode
            ? context.identity_transform
            : makeTransformFromAlignedSequence(context.base_reference_seq);

        if (!ref_gap_pos.empty() && ref_gap_pos.size() != context.baseLength()) {
            throw std::runtime_error("Reference coordinate length mismatch between consensus and aligned references");
        }

        context.ref_transforms.reserve(ref_aligned_map.size());
        for (auto& [ref_id, ref_cigar] : ref_aligned_map) {
            context.ref_transforms.emplace(ref_id, makeTransformFromCigar(std::move(ref_cigar)));
        }

        if (!options.insertion_tsv_path.empty() &&
            options.insertion_tsv_path.lexically_normal() == output.lexically_normal()) {
            throw std::runtime_error("Insertion TSV output must be different from final FASTA output");
        }

        spdlog::info("Writing final MSA FASTA: {}", output.string());
        seq_io::SeqWriter final_writer(output, U_MAX);
        WriteState state;

        const std::vector<FilePath> insertion_paths = outs_with_insertion_path;
        const std::vector<FilePath> normal_paths = outs_path;
        const std::vector<FilePath> all_query_paths = concatPaths(insertion_paths, normal_paths);

        writeInsertionTsv(insertion_paths, context, options.insertion_tsv_path);

        const auto write_references_if_requested =
            [&](const std::vector<std::size_t>& insertion_widths) {
                if (!options.write_reference) {
                    spdlog::info("Reference output disabled; skipping reference records");
                    return;
                }
                writeReferenceRecords(final_writer, consensus_aligned_file, context,
                                      insertion_widths, state, progress);
            };

        if (options.keep_length) {
            spdlog::info("Merge mode: keep-length reference projection");
            const std::vector<std::size_t> insertion_widths = zeroInsertionWidths(context.baseLength());
            write_references_if_requested(insertion_widths);
            writeProjectedSamPaths(all_query_paths, context, insertion_widths,
                                   false, options.batch_size, threads,
                                   final_writer, state, progress);
        } else if (options.insertion_merge_mode == InsertionMergeMode::reference_guided) {
            spdlog::info("Merge mode: reference-guided per-slot insertion MSA");
            const ReferenceGuidedInsertionMsa insertion_msa =
                alignReferenceGuidedInsertionSlots(
                    insertion_paths, context, reference_guided_insertion_dir,
                    msa_tool, threads);
            spdlog::info("Reference-guided insertion MSA: {} records, {} insertion segments, {} externally aligned slots, final length {}",
                         insertion_msa.records,
                         insertion_msa.segments,
                         insertion_msa.aligned_slots,
                         expandedLength(context.baseLength(), insertion_msa.widths));

            write_references_if_requested(insertion_msa.widths);
            writeProjectedSamPathsWithAlignedInsertions(
                insertion_paths, context, insertion_msa,
                options.batch_size, threads,
                final_writer, state, progress);
            writeProjectedSamPaths(normal_paths, context, insertion_msa.widths,
                                   false, options.batch_size, threads,
                                   final_writer, state, progress);
        } else {
            spdlog::info("Merge mode: external insertion MSA");
            const std::size_t insertion_count =
                writeExternalMsaInput(insertion_paths, context, insertion_fasta_path);

            if (insertion_count == 0) {
                spdlog::info("No insertion-bearing records found; external insertion MSA skipped");
                const std::vector<std::size_t> insertion_widths = zeroInsertionWidths(context.baseLength());
                write_references_if_requested(insertion_widths);
                writeProjectedSamPaths(normal_paths, context, insertion_widths,
                                       false, options.batch_size, threads,
                                       final_writer, state, progress);
            } else {
                alignConsensusSequence(insertion_fasta_path, aligned_insertion_fasta,
                                      msa_tool, threads);
                const seq_io::SeqRecord aligned_reference =
                    readFirstFastaRecord(aligned_insertion_fasta);
                const std::vector<std::size_t> insertion_widths =
                    insertionWidthsFromAlignedReference(
                        context.base_reference_seq, aligned_reference.seq);

                write_references_if_requested(insertion_widths);
                writeAlignedInsertionRecords(final_writer, aligned_insertion_fasta,
                                             state, progress);
                writeProjectedSamPaths(normal_paths, context, insertion_widths,
                                       false, options.batch_size, threads,
                                       final_writer, state, progress);
            }
        }

        final_writer.flush();
        progress.done();
        spdlog::info("Merge completed: {} sequences total, length {}",
                    state.seq_count, state.expected_length);
    }

    void RefAligner::mergeAlignedResults(const FilePath output, std::size_t batch_size)
    {
        MergeOptions options;
        options.batch_size = batch_size;
        mergeAlignedResults(output, options);
    }

    void RefAligner::removeRefGapColumns(
        std::string& seq,
        const std::vector<bool>& ref_gap_pos)
    {
        if (ref_gap_pos.empty()) {
            return;
        }

#ifdef _DEBUG
        if (seq.size() != ref_gap_pos.size()) {
            throw std::runtime_error(
                "removeRefGapColumns: sequence length mismatch seq_len=" +
                std::to_string(seq.size()) +
                ", ref_gap_pos_len=" + std::to_string(ref_gap_pos.size()));
        }
#endif

        std::size_t write_pos = 0;
        const std::size_t n = seq.size();

        for (std::size_t read_pos = 0; read_pos < n; ++read_pos) {
            if (!ref_gap_pos[read_pos]) {
                seq[write_pos++] = seq[read_pos];
            }
        }

        seq.resize(write_pos);
    }

} // namespace align
