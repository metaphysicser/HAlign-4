#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <optional>

// 程序入口：命令行解析 -> 预处理 -> 共识对齐 -> 序列比对 -> 结果合并 -> 清理工作目录

static FilePath generatedUngappedReferencePath(const Options& opt) {
    return FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / "reference_ungapped.fasta";
}

static void writeUngappedReferenceFasta(const FilePath& aligned_fasta,
                                        const FilePath& ungapped_fasta,
                                        const std::string& source_label) {
    file_io::ensureParentDirExists(ungapped_fasta);

    seq_io::KseqReader reader(aligned_fasta);
    seq_io::SeqWriter writer(ungapped_fasta);
    seq_io::SeqRecord rec;
    std::optional<std::size_t> expected_aligned_len;
    std::size_t sequence_count = 0;

    while (reader.next(rec)) {
        ++sequence_count;

        if (!expected_aligned_len.has_value()) {
            expected_aligned_len = rec.seq.size();
        } else if (rec.seq.size() != expected_aligned_len.value()) {
            throw std::runtime_error(source_label + " is not aligned: sequence '" + rec.id +
                                     "' has length " + std::to_string(rec.seq.size()) +
                                     ", expected " + std::to_string(expected_aligned_len.value()));
        }

        seq_io::cleanSequence(rec.seq);

        std::string ungapped;
        ungapped.reserve(rec.seq.size());
        for (char ch : rec.seq) {
            if (ch != '-') {
                ungapped.push_back(ch);
            }
        }
        if (ungapped.empty()) {
            throw std::runtime_error(source_label + " contains an empty reference after removing gaps for ID '" +
                                     rec.id + "'");
        }

        rec.seq = std::move(ungapped);
        writer.write(rec);
    }

    writer.flush();

    if (sequence_count == 0) {
        throw std::runtime_error(source_label + " contains no FASTA records");
    }

    spdlog::info("Generated ungapped reference FASTA from {}: {} -> {} ({} sequences)",
                 source_label, aligned_fasta.string(), ungapped_fasta.string(), sequence_count);
}

static void normalizeAlignedReferenceOptions(Options& opt) {
    if (opt.reference_is_aligned) {
        const FilePath aligned_reference = FilePath(opt.reference_path);
        const FilePath ungapped_reference = generatedUngappedReferencePath(opt);
        writeUngappedReferenceFasta(aligned_reference, ungapped_reference,
                                    "-r/--reference with -a/--reference-aligned");
        opt.reference_path = ungapped_reference.string();
        opt.reference_msa_path = aligned_reference.string();
        return;
    }

    if (opt.reference_path.empty() && !opt.reference_msa_path.empty()) {
        const FilePath aligned_reference = FilePath(opt.reference_msa_path);
        const FilePath ungapped_reference = generatedUngappedReferencePath(opt);
        writeUngappedReferenceFasta(aligned_reference, ungapped_reference, "--reference-msa");
        opt.reference_path = ungapped_reference.string();
    }
}

// 参数校验与工作目录准备
static void checkOption(Options& opt) {
    // 文件校验
    file_io::requireRegularFile(opt.input, "input");
    if (!opt.reference_path.empty()) {
        file_io::requireRegularFile(opt.reference_path, "reference_path");
    }
    if (!opt.reference_msa_path.empty()) {
        file_io::requireRegularFile(opt.reference_msa_path, "reference_msa_path");
    }
    if (opt.reference_is_aligned && opt.reference_path.empty()) {
        throw std::runtime_error("-a/--reference-aligned requires -r/--reference");
    }
    if (opt.reference_is_aligned && !opt.reference_msa_path.empty()) {
        throw std::runtime_error("-a/--reference-aligned cannot be used together with --reference-msa; "
                                 "use either '-r <aligned.fa> -a' or '-r <ungapped.fa> --reference-msa <aligned.fa>'");
    }
    if (!opt.reference_is_aligned && !opt.reference_path.empty() && !opt.reference_msa_path.empty()) {
        // 兼容旧用法：验证 -r/--reference 和 --reference-msa 的序列一致性（删除 gap 后）
        validateRefAlignedConsistency(FilePath(opt.reference_path), FilePath(opt.reference_msa_path));
    }
    if (!opt.score_matrix_path.empty()) {
        file_io::requireRegularFile(opt.score_matrix_path, "score_matrix");
        opt.score_matrix = readScoreMatrixFile(opt.score_matrix_path);
    }

    // 数值校验
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    if (opt.minimizer_size <= 0) throw std::runtime_error("minimizer-size must be > 0");
    if (opt.sketch_kmer_size <= 0) throw std::runtime_error("sketch-kmer-size must be > 0");
    if (opt.minimizer_window <= 0) throw std::runtime_error("minimizer-window must be > 0");
    if (opt.consensus_num <= 0) throw std::runtime_error("consensus-num must be > 0");
    if (opt.gap_open < 0 || opt.gap_open > 127) throw std::runtime_error("gap_open must be in [0, 127]");
    if (opt.gap_extend < 0 || opt.gap_extend > 127) throw std::runtime_error("gap_extend must be in [0, 127]");
    if (opt.min_profile_references <= 0) throw std::runtime_error("min-profile-references must be > 0");
    if (opt.max_profile_references < opt.min_profile_references) {
        throw std::runtime_error("max-profile-references must be >= min-profile-references");
    }
    if (opt.min_profile_reference_similarity < 0.0 || opt.min_profile_reference_similarity > 1.0) {
        throw std::runtime_error("min-profile-reference-similarity must be in [0, 1]");
    }
    if (!opt.insertions_output.empty() &&
        FilePath(opt.insertions_output).lexically_normal() == FilePath(opt.output).lexically_normal()) {
        throw std::runtime_error("--insertions-output must be different from -o/--output");
    }
    (void)parseInsertionMergeMode(opt.insertion_merge);
    if (opt.minimizer_size > 31) throw std::runtime_error("minimizer-size too large (must be <= 31)");
    if (opt.sketch_kmer_size > 31) throw std::runtime_error("sketch-kmer-size too large (must be <= 31)");
    if (opt.minimizer_window >= 256) {
        spdlog::warn("minimizer-window >= 256 may be slow; current value: {}", opt.minimizer_window);
    }
	// --enable-wfa 只有在 seq2seq 模式才能开启
	if (opt.enable_wfa && !opt.seq2seq)
	{
		spdlog::error("--enable-wfa can only be used with --seq2seq mode");
	}

    // workdir 准备
#ifdef _DEBUG
    constexpr bool must_be_empty = false;
#else
    constexpr bool must_be_empty = true;
#endif
    file_io::prepareEmptydir(opt.workdir, must_be_empty);

    normalizeAlignedReferenceOptions(opt);

    // MSA 命令模板解析与自检
    const std::string msa_tool_template = resolveMsaToolTemplate(opt.msa_tool);
    if (cmd::testCommandTemplate(msa_tool_template, opt.workdir, opt.threads)) {
        spdlog::info("msa-tool template test passed.");
    } else {
        throw std::runtime_error("msa-tool template test failed.");
    }
    opt.msa_tool = msa_tool_template;
}

// 清理工作目录
static void cleanupWorkdir(const Options& opt) {
    if (!opt.save_workdir) {
        try {
            spdlog::info("Removing working directory: {}", opt.workdir);
            file_io::removeAll(FilePath(opt.workdir));
            spdlog::info("Working directory removed successfully");
        } catch (const std::exception& e) {
            spdlog::warn("Failed to remove working directory: {}", e.what());
        }
    } else {
        spdlog::info("Keeping working directory: {}", opt.workdir);
    }
}

static std::size_t inferAlignmentBatchSize(uint_t sequence_count) {
    constexpr std::size_t min_batch_size = 1000;
    constexpr std::size_t max_batch_size = 10000;
    const std::size_t estimated = static_cast<std::size_t>(sequence_count) / 64U;
    return std::clamp(estimated, min_batch_size, max_batch_size);
}

int main(int argc, char** argv) {
    try
    {
        // 初始化日志
        spdlog::init_thread_pool(8192, 1);
        setupLogger();

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);
        app.formatter(std::make_shared<CustomFormatter>());
        CLI11_PARSE(app, argc, argv);

        if (opt.enable_wfa) {
            opt.seq2seq = true;
        }

        // 设置默认工作目录
        if (opt.workdir.empty()) {
            opt.workdir = makeDefaultWorkdir();
            spdlog::info("--workdir not provided, using default: {}", opt.workdir);
        }

        // 打印参数
        logParsedOptions(opt);
        spdlog::info("Starting halign4 version {}...", VERSION);

        // 校验参数
        checkOption(opt);
        setupLoggerWithFile(opt.workdir);

        // 预处理
        const uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.consensus_num);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // 文件路径定义
        const FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED;
        const FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
        const FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

        // 处理参考序列
        if (!opt.reference_path.empty())
        {
            spdlog::info("Using user-specified reference sequence: {}", opt.reference_path);
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }
            file_io::copyFile(FilePath(opt.reference_path), consensus_unaligned_file);
            spdlog::info("Reference sequence copied to: {}", consensus_unaligned_file.string());
        }

        // 快速路径：序列数 <= consensus_num 且不保留长度时直接输出。
        // 需要过滤参考或输出插入 TSV 时，仍走 RefAligner 合并路径，保证输出选项生效。
        if (preproc_count <= opt.consensus_num && opt.keep_length == false &&
            !opt.no_reference_output && opt.insertions_output.empty())
        {
            if (!opt.reference_msa_path.empty()) {
                spdlog::info("Using pre-aligned reference MSA directly: {}", opt.reference_msa_path);
                file_io::copyFile(FilePath(opt.reference_msa_path), FilePath(opt.output));
            } else {
                alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_tool, opt.threads);
                file_io::copyFile(consensus_aligned_file, FilePath(opt.output));
            }
            spdlog::info("All sequences processed; final output written to {}", opt.output);

            cleanupWorkdir(opt);
            spdlog::info("halign4 End!");
            return 0;
        }
        else if (opt.reference_path.empty())
        {
            // 生成共识序列
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_tool, opt.threads);

            const std::string consensus_string = consensus::generateConsensusSequence(
                consensus_aligned_file,
                consensus_file,
                consensus_json_file,
                opt.consensus_num,
                opt.threads,
                4096
            );

            spdlog::info("Consensus sequence generated with length {}", consensus_string.size());
        }

        // 比对阶段
        const FilePath reference_path = opt.reference_path.empty() ? consensus_file : FilePath(opt.reference_path);
        align::RefAligner ref_aligner(opt, reference_path);

        // 批大小策略：
        // - 用户显式传 --batch-size 时使用用户值；
        // - 未传时按输入序列数估计，并限制在 [1000, 10000]。
        const std::size_t cli_batch_size = (opt.batch_size > 0)
            ? static_cast<std::size_t>(opt.batch_size)
            : 0U;
        const std::size_t inferred_batch_size = inferAlignmentBatchSize(preproc_count);
        const std::size_t alignment_batch_size =
            (cli_batch_size > 0) ? cli_batch_size : inferred_batch_size;
        const std::size_t merge_batch_size = alignment_batch_size;

        // 默认走 seq2profile；仅当用户显式开启 --seq2seq 时切换到 seq2seq。
        if (opt.seq2seq) {
            spdlog::info("Alignment mode: seq2seq, batch_size={}", alignment_batch_size);
            ref_aligner.alignSeq2Seq(opt.input, alignment_batch_size);
        } else {
            spdlog::info("Alignment mode: seq2profile, batch_size={}", alignment_batch_size);
            ref_aligner.alignSeq2Profile(opt.input, alignment_batch_size);
        }

        align::MergeOptions merge_options;
        merge_options.batch_size = merge_batch_size;
        merge_options.keep_length = opt.keep_length;
        merge_options.write_reference = !opt.no_reference_output;
        merge_options.insertion_tsv_path = FilePath(opt.insertions_output);
        merge_options.insertion_merge_mode = parseInsertionMergeMode(opt.insertion_merge);
        ref_aligner.mergeAlignedResults(opt.output, merge_options);

        cleanupWorkdir(opt);

        spdlog::info("halign4 End!");
        return 0;
    } catch (const std::exception &e) {
        spdlog::error("Fatal error: {}", e.what());
        spdlog::error("halign4 End!");
        return 1;
    } catch (...) {
        spdlog::error("Fatal error: unknown exception");
        spdlog::error("halign4 End!");
        return 1;
    }
}
