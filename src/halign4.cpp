#include <config.hpp>
#include <utils.h>
#include "preprocess.h"
#include "consensus.h"

#include "align.h"

// 该文件是程序入口（main），负责：
// 1) 解析命令行参数（使用 CLI11），并将值绑定到 Options 结构；
// 2) 校验参数（例如文件存在性、工作目录准备、外部 MSA 命令模板测试）；
// 3) 调用预处理（preprocessInputFasta）将原始输入复制/下载到工作目录并清洗；
// 4) 选择/生成参考序列（center_path 或共识序列）并完成比对与结果合并；
// 5) 输出最终的对齐文件，并按需清理工作目录。
//
// 注：本文件尽量把流程拆成职责清晰的小步骤，便于定位错误并在 CI 中逐步测试。

// checkOption：参数再校验与工作目录准备
// 该函数做了多项关键工作：
// - 使用 file_io 模块统一校验文件存在性与路径类型（统一错误信息与行为）；
// - 检查数值参数的合理性（threads/kmer_size/kmer_window/cons_n）；
// - 准备工作目录：Debug 模式允许非空 workdir 以便快速迭代；Release 模式要求 workdir 为空以避免覆盖已有数据；
// - 对 MSA 命令模板做一次自检（cmd::testCommandTemplate）：在临时目录运行模板并验证是否产生输出，用于尽早发现模板占位符/命令行问题。
//   注意：testCommandTemplate 会实际运行外部命令，可能依赖外部工具或运行环境；若不希望在参数校验阶段执行外部命令，可将测试延后并在错误时给出更明确提示。
static void checkOption(Options& opt) {
    // 文件：输入必须存在且为普通文件
    file_io::requireRegularFile(opt.input, "input");

    // 文件：若指定 center_path，则也必须存在
    if (!opt.center_path.empty()) {
        file_io::requireRegularFile(opt.center_path, "center_path");
    }

    // 数值：线程数必须为正
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    // 数值：k 必须为正
    if (opt.kmer_size <= 0) throw std::runtime_error("kmer_size must be > 0");
    // 数值：窗口大小必须为正
    if (opt.kmer_window <= 0) throw std::runtime_error("kmer_window must be > 0");
    // 数值：共识使用的 Top-N 必须为正
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");

    // k 的上限约束：避免超出下游实现的设计范围（保持现有逻辑）
    if (opt.kmer_size > 31) throw std::runtime_error("kmer_size too large (must be <= 31)");
    // w 过大可能导致性能下降或走非优化路径：这里仅告警，不改变行为
    if (opt.kmer_window >= 256) {
        spdlog::warn("kmer_window >= 256 may be slow/unsupported by the fast minimizer path; current value: {}", opt.kmer_window);
    }

    // workdir：Release 下必须为空目录（防止覆盖）；Debug 下允许复用非空目录便于调试
#ifdef _DEBUG
    constexpr bool must_be_empty = false;
#else
    constexpr bool must_be_empty = true;
#endif
    // 创建/校验工作目录（必要时清空；具体策略由 must_be_empty 控制）
    file_io::prepareEmptydir(opt.workdir, must_be_empty);

    // msa_cmd：允许输入关键字（minipoa/mafft/clustalo）或自定义“命令模板字符串”
    // 这里将关键字解析成最终模板，避免后续模块重复分支判断
    const std::string msa_cmd_str = resolveMsaCmdTemplate(opt.msa_cmd);

    // 模板自检：在工作目录中跑一个小样例，尽早暴露外部命令不可用/占位符错误等问题
    if (cmd::testCommandTemplate(msa_cmd_str, opt.workdir, opt.threads)) {
        spdlog::info("msa_cmd template test passed.");
    } else {
        // 自检失败直接终止：避免后续跑到中途才报错（保持现有行为）
        throw std::runtime_error("msa_cmd template test failed.");
    }

    // 规范化：后续模块统一使用“最终模板字符串”
    opt.msa_cmd = msa_cmd_str;

    // 可选：确保输出父目录存在（当前保持现有逻辑：不自动创建）
    // file_io::ensureParentDirExists(opt.output);
}

// cleanupWorkdir：根据 save_workdir 开关决定是否删除工作目录
// 说明：
// - 当 save_workdir=false（默认）时，程序在成功完成后会尝试删除 workdir，以避免堆积中间文件；
// - 当 save_workdir=true 时，保留 workdir 便于 debug/复现/检查中间产物；
// - 删除失败不会影响主流程返回值：这里用 warn 记录原因，避免“成功结果被清理失败掩盖”。
static void cleanupWorkdir(const Options& opt) {
    // 开关：默认删除，除非用户显式要求保留
    if (!opt.save_workdir) {
        try {
            // 日志：记录即将删除的目录
            spdlog::info("Removing working directory: {}", opt.workdir);
            // 递归删除：清理所有中间产物
            file_io::removeAll(FilePath(opt.workdir));
            // 日志：删除成功
            spdlog::info("Working directory removed successfully");
        } catch (const std::exception& e) {
            // 清理失败不影响主流程：仅告警
            spdlog::warn("Failed to remove working directory: {}", e.what());
        }
    } else {
        // 保留工作目录：便于定位问题
        spdlog::info("Keeping working directory: {}", opt.workdir);
    }
}

int main(int argc, char** argv) {
    // main 的整体流程（高层概览）：
    // 1) 初始化日志与线程池
    // 2) 定义 CLI 并解析参数
    // 3) 打印并校验参数
    // 4) 将日志输出到工作目录
    // 5) 执行预处理（fetch/copy + 清洗 + Top-K 选择）
    // 6) 决定参考序列来源：用户 center_path 或预处理生成的共识序列
    // 7) 快速路径：当序列数 <= cons_n 且不要求保留长度策略时，直接 MSA 并输出
    // 8) 否则使用 RefAligner 比对（alignQueryToRef）并合并（mergeAlignedResults）
    // 9) 输出最终对齐结果
    // 10) 根据 --save-workdir 决定是否清理 workdir

    try
    {
        // spdlog：初始化异步日志线程池（队列大小 8192，后端线程 1）
        spdlog::init_thread_pool(8192, 1);
        // 默认日志配置（控制台等）
        setupLogger();

        // 保存命令行参数的结构体
        Options opt;
        // CLI11：创建命令行解析器
        CLI::App app{"halign4"};

        // 注册命令行参数到 opt（定义 --input/--output/--threads 等）
        setupCli(app, opt);

        // CLI11：解析 argc/argv（失败将自动打印帮助并退出）
        CLI11_PARSE(app, argc, argv);

        // 若用户未提供 -w/--workdir，则在此处生成默认工作目录。
        // 这样做能保证：
        // 1) CLI 帮助信息不会出现“每次不同的随机默认值”；
        // 2) checkOption 中仍可统一调用 file_io::prepareEmptydir 创建/校验目录；
        // 3) 程序后续所有路径拼接逻辑保持不变（只是 workdir 有了默认值）。
        if (opt.workdir.empty()) {
            // 生成默认 workdir（通常含时间戳/随机后缀，避免冲突）
            opt.workdir = makeDefaultWorkdir();
            // 记录默认 workdir 供用户排查
            spdlog::info("--workdir not provided, using default workdir: {}", opt.workdir);
        }

        // 打印解析后的参数（便于复现与排错）
        logParsedOptions(opt);
        // 程序启动：打印版本号
        spdlog::info("Starting halign4 version {}...", VERSION);

        // 参数校验（文件存在性、workdir、msa 模板自检等）
        checkOption(opt);

        // 将日志输出到 workdir（便于收集运行时日志）
        // 注意：该调用依赖 workdir 已在 checkOption 中准备好
        setupLoggerWithFile(opt.workdir);

        // ---------------- 预处理阶段（IO 密集） ----------------
        // preprocessInputFasta 的职责：
        // - 如果 input 是 URL 则下载到 workdir/data/raw，否则拷贝到该目录；
        // - 逐条读取、清洗并写出到 data/clean，同时选出最长的 cons_n 条序列写入共识输入文件；
        // - 返回处理的总序列数（用于决定下一步是否走快速路径/是否需要 RefAligner 合并）。
        const uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.cons_n);
        // 日志：记录总条数
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // ---------------- 共识对齐/参考准备阶段（可能调用外部 MSA） ----------------
        // 工作目录中关键中间文件的路径约定（由 preprocess/consensus 模块负责写入/读取）
        const FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_UNALIGNED; // 未对齐的共识输入
        const FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_ALIGNED;     // 外部 MSA 的对齐输出
        const FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_FASTA;               // 生成的共识序列（FASTA）
        const FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN / CLEAN_CONS_JSON;           // 共识统计/元信息（JSON）

        if (!opt.center_path.empty())
        {
            // 用户指定 center_path：将其拷贝为 workdir 内的共识输入文件，保持后续流程统一
            // 说明：
            // 1) 不修改用户原始文件；
            // 2) 若 workdir 内已有旧的 consensus_unaligned_file，则先删除再拷贝，避免混用旧数据。
            spdlog::info("Using user-specified center sequence: {}", opt.center_path);

            // 删除工作目录中旧文件（如果存在）
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }

            // 拷贝 center_path -> workdir 的 consensus_unaligned_file
            file_io::copyFile(FilePath(opt.center_path), consensus_unaligned_file);
            spdlog::info("Center sequence copied to: {}", consensus_unaligned_file.string());
        }

        // 快速路径：如果总序列数 <= cons_n，且不启用 keep_length，
        // 则无需 RefAligner 的“参考比对 + 合并”，直接对 cons_n 子集做一次 MSA 即可输出。
        if (preproc_count <= opt.cons_n && opt.keep_length == false)
        {
            // 调用外部 MSA：consensus_unaligned_file -> consensus_aligned_file
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);
            // 最终输出：直接拷贝对齐结果到用户指定 output
            file_io::copyFile(consensus_aligned_file, FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);

            // 成功退出前按需清理 workdir
            cleanupWorkdir(opt);

            spdlog::info("halign4 End!");
            return 0;
        }
        else if (opt.center_path.empty())
        {
            // 非快速路径且用户未指定 center_path：需要先对共识子集做 MSA，再从对齐结果生成共识序列。
            // TODO：未来可增加参数允许用户直接提供“已对齐的共识文件”，跳过外部 MSA。

            // 外部 MSA：得到用于生成共识的多序列对齐
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.threads);

            // 从对齐后的共识子集生成“单条参考序列”（写入 consensus_file，并产出 JSON 统计）
            const std::string consensus_string = consensus::generateConsensusSequence(
                consensus_aligned_file, // 输入：对齐后的共识子集
                consensus_file,         // 输出：共识序列 FASTA
                consensus_json_file,    // 输出：共识统计 JSON
                opt.cons_n,             // 参数：使用的序列条数
                opt.threads,            // 参数：线程数（内部可能 OpenMP 并行）
                4096                    // 参数：batch_size（读写批大小）
            );

            // 日志：打印共识序列长度（便于观察是否异常过短/过长）
            spdlog::info("Consensus sequence generated with length {}", consensus_string.size());
        }

        // ---------------- 比对阶段（使用 RefAligner） ----------------
        // 参考序列路径选择：
        // - 若 center_path 为空，则使用刚生成的 consensus_file；
        // - 否则使用用户传入的 center_path（注意：此处不强制使用 workdir 内的副本，保持现有逻辑不变）。
        const FilePath ref_path = opt.center_path.empty() ? consensus_file : FilePath(opt.center_path);

        // 创建对齐器：内部加载参考序列并准备索引
        align::RefAligner ref_aligner(opt, ref_path);
        // 比对：将输入查询序列对齐到参考序列
        ref_aligner.alignQueryToRef(opt.input);
        // 合并：汇总各块结果并按需要再次调用 MSA 做增量合并/修整，然后写出最终 output
        ref_aligner.mergeAlignedResults(opt.output, 25600);

        // 成功完成后，根据 --save-workdir 决定是否清理 workdir
        cleanupWorkdir(opt);

        // 正常结束
        spdlog::info("halign4 End!");

        return 0;
    } catch (const std::exception &e) {
        // 捕获并打印标准异常
        spdlog::error("Fatal error: {}", e.what());
        spdlog::error("halign4 End!");
        return 1;
    } catch (...) {
        // 捕获未知异常
        spdlog::error("Fatal error: unknown exception");
        spdlog::error("halign4 End!");
        return 1;
    }
}
