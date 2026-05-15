#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"

#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <memory>

// 程序入口：命令行解析 -> 预处理 -> 共识对齐 -> 序列比对 -> 结果合并 -> 清理工作目录

// 参数校验与工作目录准备
static void checkOption(Options& opt) {
    // 文件校验
    file_io::requireRegularFile(opt.input, "input");
    if (!opt.ref_path.empty()) {
        file_io::requireRegularFile(opt.ref_path, "ref_path");
    }
    if (!opt.ref_align_path.empty()) {
        file_io::requireRegularFile(opt.ref_align_path, "ref_align_path");
        if (opt.ref_path.empty()) {
            throw std::runtime_error("--ref-align requires -r/--ref to be provided as the matching reference FASTA");
        }
        // 验证 -r 和 --ref-align 的序列一致性（删除 gap 后）
        validateRefAlignedConsistency(FilePath(opt.ref_path), FilePath(opt.ref_align_path));
    }
    if (!opt.score_path.empty()) {
        file_io::requireRegularFile(opt.score_path, "score");
        opt.score_matrix = readScoreMatrixFile(opt.score_path);
    }

    // 数值校验
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    if (opt.kmer_size <= 0) throw std::runtime_error("kmer_size must be > 0");
    if (opt.profile_ref_kmer_len <= 0) throw std::runtime_error("profile_ref_kmer_len must be > 0");
    if (opt.kmer_window <= 0) throw std::runtime_error("kmer_window must be > 0");
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");
    if (opt.gap_open < 0 || opt.gap_open > 127) throw std::runtime_error("gap_open must be in [0, 127]");
    if (opt.gap_extend < 0 || opt.gap_extend > 127) throw std::runtime_error("gap_extend must be in [0, 127]");
    if (opt.profile_ref_min <= 0) throw std::runtime_error("profile_ref_min must be > 0");
    if (opt.profile_ref_max < opt.profile_ref_min) {
        throw std::runtime_error("profile_ref_max must be >= profile_ref_min");
    }
    if (opt.profile_ref_min_similarity < 0.0 || opt.profile_ref_min_similarity > 1.0) {
        throw std::runtime_error("profile_ref_min_similarity must be in [0, 1]");
    }
    if (opt.kmer_size > 31) throw std::runtime_error("kmer_size too large (must be <= 31)");
    if (opt.profile_ref_kmer_len > 31) throw std::runtime_error("profile_ref_kmer_len too large (must be <= 31)");
    if (opt.kmer_window >= 256) {
        spdlog::warn("kmer_window >= 256 may be slow; current value: {}", opt.kmer_window);
    }
	// --wfa只有在seq2seq模式才能开启
	if (opt.wfa && !opt.seq2seq)
	{
		spdlog::error("--wfa can only be used with --seq2seq mode");
	}

    // workdir 准备
#ifdef _DEBUG
    constexpr bool must_be_empty = false;
#else
    constexpr bool must_be_empty = true;
#endif
    file_io::prepareEmptydir(opt.workdir, must_be_empty);

    // MSA 命令模板解析与自检
    const std::string msa_cmd_str = resolveMsaCmdTemplate(opt.msa_cmd);
    if (cmd::testCommandTemplate(msa_cmd_str, opt.workdir, opt.threads)) {
        spdlog::info("msa_cmd template test passed.");
    } else {
        throw std::runtime_error("msa_cmd template test failed.");
    }
    opt.msa_cmd = msa_cmd_str;
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
        const uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.cons_n);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // 文件路径定义
        const FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED;
        const FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;
        const FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;
        const FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;

        // 处理参考序列
        if (!opt.ref_path.empty())
        {
            spdlog::info("Using user-specified reference sequence: {}", opt.ref_path);
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }
            file_io::copyFile(FilePath(opt.ref_path), consensus_unaligned_file);
            spdlog::info("Reference sequence copied to: {}", consensus_unaligned_file.string());
        }

        // 快速路径：序列数 <= cons_n 且不保留长度时直接输出
        if (preproc_count <= opt.cons_n && opt.keep_length == false)
        {
            if (!opt.ref_align_path.empty()) {
                spdlog::info("Using pre-aligned reference MSA directly: {}", opt.ref_align_path);
                file_io::copyFile(FilePath(opt.ref_align_path), FilePath(opt.output));
            } else {
                alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);
                file_io::copyFile(consensus_aligned_file, FilePath(opt.output));
            }
            spdlog::info("All sequences processed; final output written to {}", opt.output);

            cleanupWorkdir(opt);
            spdlog::info("halign4 End!");
            return 0;
        }
        else if (opt.ref_path.empty())
        {
            // 生成共识序列
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);

            const std::string consensus_string = consensus::generateConsensusSequence(
                consensus_aligned_file,
                consensus_file,
                consensus_json_file,
                opt.cons_n,
                opt.threads,
                4096
            );

            spdlog::info("Consensus sequence generated with length {}", consensus_string.size());
        }

        // 比对阶段
        const FilePath ref_path = opt.ref_path.empty() ? consensus_file : FilePath(opt.ref_path);
        align::RefAligner ref_aligner(opt, ref_path);

        // 批大小策略：
        // - 用户显式传 --batch-size 时，统一覆盖 align 与 merge 的 batch；
        // - 未传（0）时保持历史行为：seq2profile=128，seq2seq=函数默认值(25600)，merge=25600。
        const std::size_t cli_batch_size = (opt.batch_size > 0)
            ? static_cast<std::size_t>(opt.batch_size)
            : 0U;
        const std::size_t merge_batch_size = (cli_batch_size > 0) ? cli_batch_size : 25600U;

        // 默认走 seq2profile；仅当用户显式开启 --seq2seq 时切换到 seq2seq。
        if (opt.seq2seq) {
            spdlog::info("Alignment mode: seq2seq, batch_size={}",
                         (cli_batch_size > 0 ? cli_batch_size : 25600U));
            ref_aligner.alignSeq2Seq(opt.input, cli_batch_size);
        } else {
            const std::size_t seq2profile_batch_size = (cli_batch_size > 0) ? cli_batch_size : 1000;
            spdlog::info("Alignment mode: seq2profile, batch_size={}", seq2profile_batch_size);
            ref_aligner.alignSeq2Profile(opt.input, seq2profile_batch_size);
        }
        ref_aligner.mergeAlignedResults(opt.output, merge_batch_size);

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
