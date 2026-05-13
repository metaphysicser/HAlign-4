#ifndef CONFIG_HPP
#define CONFIG_HPP

// ------------------------------------------------------------------
// config.hpp
// 说明（详细中文注释）
//
// 本文件集中定义了全局配置常量、日志初始化函数、以及与命令行解析（CLI11）相关的
// 辅助类型与工具（例如自定义格式器与 validator）。
//
// 目的：
// - 为项目提供单一入口的“配置库”，便于在代码各处引用一致的符号（例如工作目录结构、默认命令模板、日志文件名等）；
// - 提供便捷的日志初始化函数 `setupLogger` / `setupLoggerWithFile`，方便在 main 中统一配置日志输出到控制台与文件；
// - 提供 CLI 美化与输入修剪（trim_whitespace），提高命令行体验与容错性。
//
// 注意事项：
// - 这里的默认 MSA 命令模板（DEFALT_MSA_CMD）仅做占位和演示用途；生产环境请替换为实际可用的 MSA 工具与参数。
// - 所有路径常量为字符串字面量（相对路径），在使用时通常与 `workdir` 进行拼接以形成绝对或工作相对路径。
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// 引入头文件：功能模块包括线程池、命令行解析、日志系统、序列化库等
// ------------------------------------------------------------------
#include <CLI/CLI.hpp>                       // CLI11 命令行解析库（本地版本）
#include "spdlog/spdlog.h"                       // spdlog 主头文件
#include "spdlog/sinks/stdout_color_sinks.h"     // 控制台彩色输出 sink
#include "spdlog/sinks/basic_file_sink.h"        // 文件输出 sink
#include "spdlog/async.h"                        // 异步日志支持

#include <filesystem>
#include <array>
#include <sstream>
#include <cinttypes>
#include <random>
#include <chrono>
#include <iomanip>

// ------------------------------------------------------------------
// 通用配置常量
// ------------------------------------------------------------------
#define VERSION "2.0.0"                   // 版本号，程序启动时可打印以便追踪
#define LOGGER_NAME "logger"              // 默认日志器名称（用于 spdlog 注册）
#define LOGGER_FILE "halign4.log"         // 默认日志文件名（相对于工作目录）
#define CONFIG_FILE "config.json"         // 默认配置文件路径（如果将来支持外部配置）

const std::string MINIPOA_CMD = "minipoa {input} -S -t {thread} -r1 > {output}"; // Minipoa 多序列比对命令模板示例
const std::string MAFFT_MSA_CMD = "mafft --thread {thread} --auto {input} > {output}"; // MAFFT 多序列比对命令模板示例
const std::string CLUSTALO_MSA_CMD = "clustalo -i {input} -o {output} --threads {thread}"; // Clustal Omega 多序列比对命令模板示例

const std::string DEFAULT_MSA_CMD = MINIPOA_CMD; // 默认多序列比对命令模板

// 默认 DNA5 打分矩阵（A/C/G/T/N），矩阵文件与内部顺序均使用该顺序。
static constexpr std::array<int8_t, 25> DEFAULT_DNA5_SCORE_MATRIX = {
     4, -2,  1, -2,  0,
    -2,  4, -2,  1,  0,
     1, -2,  4, -2,  0,
    -2,  1, -2,  4,  0,
     0,  0,  0,  0,  0
};

// ------------------------------------------------------------------
// resolveMsaCmdTemplate：把用户在 -p/--msa-cmd 中输入的内容解析成“最终命令模板”。
//
// 需求：
// - 当用户输入 minipoa / mafft / clustalo 时，自动使用对应的内置模板命令；
// - 当用户输入的是自定义模板（包含 {input}/{output} 等占位符）时，保持原样；
// - 当用户不输入 -p 时，沿用 DEFAULT_MSA_CMD，不改变现有默认行为。

// 设计说明（正确性/可用性）：
// - 不能再把 -p 当作“文件路径”去校验（ExistingFile / requireRegularFile），因为这些工具名通常依赖 PATH。
// - 这里只对“完全等于关键字”的情况做映射，避免误伤用户自定义命令（例如 "mafft --auto ..."）。
// ------------------------------------------------------------------
inline std::string resolveMsaCmdTemplate(const std::string& user_value) {
    // trim：去除两端空白，避免用户误输入空格导致关键字匹配失败
    const auto start = user_value.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) {
        return DEFAULT_MSA_CMD;
    }
    const auto end = user_value.find_last_not_of(" \t\n\r");
    std::string v = user_value.substr(start, end - start + 1);

    // tolower：关键字大小写不敏感
    for (char& c : v) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (v == "minipoa") return MINIPOA_CMD;
    if (v == "mafft") return MAFFT_MSA_CMD;
    if (v == "clustalo") return CLUSTALO_MSA_CMD;

    // 其他情况：认为是用户自定义模板，原样返回
    return user_value;
}

// 工作目录体系
const std::string WORKDIR_DATA = "data";         // 原始数据目录
const std::string WORKDIR_TMP = "temp";         // 临时目录（短生命周期文件）
const std::string RESULTS_DIR = "result";      // 最终结果目录

const std::string DATA_RAW = "raw_data";        // 原始数据子目录（下载/拷贝来的原始输入）
const std::string DATA_CLEAN = "clean_data";    // 清理后数据子目录（预处理后写出的序列）

// 共识相关的文件名（相对于 DATA_CLEAN）
const std::string CLEAN_CONS_UNALIGNED = "consensus_unaligned.fasta"; // 共识序列文件名（未对齐）
const std::string CLEAN_CONS_ALIGNED = "consensus_aligned.fasta";     // 共识序列文件名（已对齐）

const std::string CLEAN_CONS_FASTA = "consensus.fasta"; // 最终共识序列 FASTA 文件名
const std::string CLEAN_CONS_JSON = "consensus.json";   // 共识统计/计数输出（JSON）

// ------------------------------------------------------------------
// 对齐输出相关的文件名（相对于 RESULTS_DIR）
// ------------------------------------------------------------------
// 说明：这些文件名用于多序列比对（MSA）和插入序列处理的中间文件与最终输出文件
// 目的：集中管理文件名常量，避免在代码中硬编码字符串字面量，便于统一修改与维护

// 最终对齐结果文件名（所有序列对齐后的 MSA 输出）
#define FINAL_ALIGNED_FASTA "final_aligned.fasta"

// 插入序列相关文件名
#define ALL_INSERTION_FASTA "all_insertion.fasta"           // 合并的所有插入序列（未对齐）
#define ALIGNED_INSERTION_FASTA "aligned_insertion.fasta"   // 对齐后的插入序列（MSA 结果）

// 线程级别输出文件名模板（用于并行写出）
// 说明：每个线程独立写出 SAM 文件，避免线程竞争；最后由主线程合并
#define THREAD_SAM_PREFIX "thread"                          // 线程 SAM 文件前缀（"thread" + tid + ".sam"）
#define THREAD_SAM_SUFFIX ".sam"                            // 线程 SAM 文件后缀
#define THREAD_INSERTION_SAM_SUFFIX "_insertion.sam"        // 线程插入序列 SAM 文件后缀（"thread" + tid + "_insertion.sam"）

// ------------------------------------------------------------------
// 调试与整数精度配置
// ------------------------------------------------------------------
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef M64
#define M64 0
#endif

// 根据宏定义切换 32 位或 64 位整数
// 目的：在处理非常大量的数据计数时，使用 64-bit 能避免溢出；开发/轻量运行可用 32-bit 节省内存。
#if M64
typedef int64_t	int_t;
typedef uint64_t uint_t;
#define PRIdN	PRId64
#define U_MAX	UINT64_MAX
#define I_MAX	INT64_MAX
#define I_MIN	INT64_MIN
#else
typedef int32_t int_t;
typedef uint32_t uint_t;
#define PRIdN	PRId32
#define U_MAX	UINT32_MAX
#define I_MAX	INT32_MAX
#define I_MIN	INT32_MIN
#endif

// 获取硬件并发线程数（兜底 1）
static int get_default_threads() {
    unsigned int hc = std::thread::hardware_concurrency();
    return static_cast<int>(hc ? hc : 1u);
}

// ------------------------------------------------------------------
// 默认 workdir 生成器（关键逻辑新增，需中文注释）
//
// 需求：用户不传 -w/--workdir 时，自动使用 "./tmp-随机数"。
// 设计点：
// 1) 这里返回的是相对路径字符串（"./tmp-..."），保持与 CLI 输入一致；
// 2) 使用“时间戳 + 随机数”拼接，降低并发/重复运行时的碰撞概率；
// 3) 不在此处创建目录：目录的创建/清空策略仍由 checkOption()->file_io::prepareEmptydir 统一处理，
//    以保证现有流程与错误处理逻辑不变。
// ------------------------------------------------------------------
static std::string makeDefaultWorkdir() {
    using Clock = std::chrono::high_resolution_clock;
    const auto now = Clock::now().time_since_epoch();
    const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();

    // 使用 random_device 作为种子，生成少量随机扰动，进一步降低同一纳秒内启动时的冲突概率。
    std::random_device rd;
    std::mt19937_64 gen(static_cast<uint64_t>(rd()) ^ static_cast<uint64_t>(ns));
    std::uniform_int_distribution<uint32_t> dist(0u, 0xFFFFFFFFu);
    const uint32_t r = dist(gen);

    std::ostringstream oss;
    oss << "./tmp-" << ns << "-" << std::hex << std::setw(8) << std::setfill('0') << r;
    return oss.str();
}

struct Options {
	// 输入/输出与工作目录
	std::string input;          // -i：输入序列文件（路径或拷贝文件）
	std::string output;         // -o：最终输出文件（写入位置）
	std::string workdir;        // -w：工作目录，所有中间文件（data/raw, data/clean 等）放在该目录下

	// 可选参数：参考序列、MSA 命令模板
	std::string ref_path;       // -r：可选，指定参考/中心序列文件路径，若指定则绕过自动选择
    std::string ref_align_path; // --ref-align：与 -r/--ref 对应的“已比对 MSA”文件路径
	std::string msa_cmd;        // -p：用于对共识序列做 MSA 的命令模板（可以包含 {input} {output} {thread} 占位符）

	// 并行与算法参数
	int threads = get_default_threads(); // -t：线程数，默认为 CPU 核心数
	int kmer_size = 19;         // --kmer-size：用于归类/聚类的 k-mer 大小（后续步骤使用）
	int kmer_window = 19;       // --kmer-window：minimizer 窗口大小 w（以 k-mer 为单位）
	int cons_n = 1000;          // --cons-n：挑选用于共识计算的序列数量（Top-K by length）
	int sketch_size = 3000;     // --sketch-size：用于 sketch 的大小（默认 3000）
    // 说明：将 sketch 的 k-mer 大小与 minimizer 的 k-mer 大小解耦。
    // - kmer_size 仍用于 minimizer/锚点；
    // - sketch_kmer_size 仅用于 mash::sketchFromSequence，默认 21。
    int sketch_kmer_size = 21;  // --sketch-kmer-size：用于 sketch 的 k-mer 大小（默认 21）
	int batch_size = 0;         // --batch-size：对齐批大小；0 表示按模式使用内置默认值
	bool wfa = false;           // --wfa：启用 WFA 开关（默认关闭，保证现有行为不变）
	bool seq2seq = false;       // --seq2seq：开启后使用 seq2seq 比对路径；默认使用 seq2profile
    std::string score_path;     // -s/--score：可选 DNA5 打分矩阵文件路径
    std::array<int8_t, 25> score_matrix = DEFAULT_DNA5_SCORE_MATRIX;
    int gap_open = 10;          // --gap-open：gap open 罚分
    int gap_extend = 2;         // --gap-extend：gap extend 罚分

	// keep length 相关开关：
	// - keep_length：保持“第一条/中心序列”的长度不变（其余序列允许按对齐结果变化/填充），适用于只关心输出共识/中心序列长度的场景。
	bool keep_length = false; // --keep-length
	// workdir 清理开关：
	// - save_workdir：若为 true，则在完成比对/输出后保留工作目录；若为 false（默认），则在成功完成后删除工作目录。
	bool save_workdir = false;      // --save-workdir
};


// setupCli：定义 CLI 参数并绑定到 Options
// 注释说明：
// - 使用 CLI11 库实现参数解析，支持短参数/长参数与基本校验（例如 ExistingFile）
// - 对于像 -p/--msa-cmd 这类可能是可执行名（而非完整路径）的参数，ExistingFile 会拒绝仅命令名的情况；
//   如果希望允许命令名（在 PATH 中解析），可以去掉 check(CLI::ExistingFile) 或改为用户层面的更宽容判断。
static void setupCli(CLI::App& app, Options& opt) {
    app.description("HAlign 4: A New Strategy for Rapidly Aligning Millions of Sequences");

    // 设置版本标志：--version 打印版本并退出，-h/--help 也会显示版本信息
    app.set_version_flag("-v,--version", std::string("halign4 version ") + VERSION);

   // 必须参数（同时支持短参数和长参数）
    // -i/--input：输入序列（FASTA）。
    // 说明：
    // - 必须是本地文件路径（当前实现会在参数校验阶段检查是否存在）；
    // - 文件格式建议为 .fasta/.fa（内部读取使用 kseq）。
    // 示例：-i test/data/mt1x/mt1x.fasta
    app.add_option("-i,--input", opt.input,
                   "Input sequences in FASTA format (local file path).")
        ->required()
        ->check(CLI::ExistingFile);

    // -o/--output：最终输出对齐结果（FASTA）。
    // 说明：
    // - 输出路径不要求预先存在；父目录建议存在（若不存在，后续写出阶段可能失败）。
    // - 输出内容为多序列对齐后的 FASTA。
    // 示例：-o out/aligned.fasta
    app.add_option("-o,--output", opt.output,
                   "Output aligned sequences (FASTA file path).")
        ->required();

    // workdir：改为可选。
    // 若用户不提供 -w，则在 main() 解析完成后用 makeDefaultWorkdir() 生成默认值。
    // 这样做的原因：CLI11 的 default_val 会直接在 -h 中展示默认值；但这里默认值带随机数，
    // 展示出来会让帮助信息“每次不同”，也会误导用户认为必须指定固定路径。
    // -w/--workdir：工作目录（存放中间文件/日志/临时结果）。
    // 说明：
    // - 若不提供，程序会自动生成 ./tmp-<随机数>；
    // - 目录下会创建 data/raw_data、data/clean_data、temp、result 等子目录；
    // - Release 模式下要求 workdir 为空目录（避免覆盖旧结果）；Debug 模式允许复用（便于迭代）。
    app.add_option("-w,--workdir", opt.workdir,
                   "Working directory for intermediate files (default: ./tmp-<random>). ")
        ->capture_default_str();

    // 可选参数（增加长参数形式）
    // 说明：
    // - 不提供时：程序会在预处理阶段自动选择并生成共识/中心序列；
    // - 提供时：会使用该序列作为参考（并在 workdir 中进行统一管理）。
    // 典型用途：COVID 数据集里可以用 covid-ref 的第一条（武汉参考）作为 center。
    app.add_option("-r,--ref", opt.ref_path,
                   "Center/reference sequence in FASTA (optional). If not set, a consensus/center is generated.")
        ->check(CLI::ExistingFile);

    // --ref-align：与 -r/--ref 配套的“已对齐 MSA”文件。
    // 设计目标：
    // - 当用户已经为 -r 准备好了比对结果时，可以直接复用该 MSA，避免再次运行 MSA 工具；
    // - 仍然保留 -r 作为“参考/中心序列 FASTA”的入口，保持旧行为不变；
    // - 该参数不是可执行文件，也不是普通的“可选字符串”，因此直接按存在文件路径校验。
    app.add_option("-a,--align-ref", opt.ref_align_path,
                   "Pre-aligned MSA file for the reference set provided by -r/--ref.")
        ->check(CLI::ExistingFile);

    // 如果 -p 是“可执行文件路径”，ExistingFile 通常也能用；
    // 若你希望允许仅命令名（在 PATH 中），这里就不要 check
    // msa-cmd：支持关键字或自定义命令模板。
    // - 关键字：minipoa / mafft / clustalo
    // - 自定义模板：例如 "mafft --thread {thread} --auto {input} > {output}"
    // 注意：这里不能使用 ExistingFile 校验，否则关键字/命令名会被错误拒绝。
    // -p/--msa-cmd：高质量 MSA 工具（用于对共识/插入序列进行高质量对齐）。
    // 支持两种形式：
    // 1) 关键字：minipoa / mafft / clustalo
    //    - 输入关键字后，程序会自动展开为内置模板（见 MINIPOA_CMD/MAFFT_MSA_CMD/CLUSTALO_MSA_CMD）；
    // 2) 自定义“命令模板字符串”：必须至少包含 {input} 和 {output}；可选包含 {thread}。
    //    - 例如："mafft --thread {thread} --auto {input} > {output}"
    // 注意：
    // - 默认使用 minipoa（不传 -p 等价于 -p minipoa）；
    // - 该命令会在参数校验阶段用一个 tiny.fasta 做一次 smoke test，若环境缺少该工具会直接报错。
    app.add_option("-p,--msa-cmd", opt.msa_cmd,
                   "High-quality MSA method: keyword {minipoa|mafft|clustalo} or a custom command template containing {input} and {output} (optional {thread}).");

    // -t/--thread：线程数。
    // 说明：
    // - 默认值为硬件并发数（std::thread::hardware_concurrency）；
    // - 影响预处理、比对以及外部 MSA 命令中的 {thread} 替换。
    app.add_option("-t,--thread", opt.threads, "Number of threads.")
        ->default_val(get_default_threads())
        ->check(CLI::Range(1, 100000));

    // --kmer-size：k-mer 大小。
    // 说明：
    // - 用于 minimizer/哈希相关流程的参数；一般无需改动。
    // - 合法范围 [4,31]（与部分位运算/编码实现约束一致）。
    app.add_option("--kmer-size", opt.kmer_size, "K-mer size used in sketch/minimizer.")
        ->default_val(19)
        ->check(CLI::Range(4, 31));


    // --kmer-window：minimizer 窗口大小 w（单位：k-mer 数）。
    // 说明：w 越大，minimizer 更稀疏；w 越小，种子更密集但可能更慢。
    app.add_option("--kmer-window", opt.kmer_window,
                   "Minimizer window size w (in number of k-mers).")
        ->default_val(19)
        ->check(CLI::Range(1, 1000000));

    // --cons-n：用于生成共识/中心序列的 Top-N（按长度挑选）。
    // 说明：
    // - 输入序列数 <= cons_n 时，程序会直接调用外部 MSA 对全部序列做一次对齐（快速路径）；
    // - 输入序列数远大于 cons_n 时，先用 Top-N 生成共识，再进行分批/参考比对。
    app.add_option("--cons-n", opt.cons_n,
                   "Number of sequences used to build the consensus/center (Top-N by length).")
        ->default_val(1000)
        ->check(CLI::Range(1, 1000000));

    // --sketch-size：sketch（minhash）大小。
    // 说明：越大越稳健但更慢/更占内存；一般默认 3000 足够。
    app.add_option("--sketch-size", opt.sketch_size, "Sketch size (minhash count).")
        ->default_val(3000)
        ->check(CLI::Range(1, 10000000));

    // --sketch-kmer-size：仅用于 sketch 构建的 k-mer 大小。
    // 说明：
    // - 该参数不会影响 minimizer 的 k（minimizer 仍使用 --kmer-size）；
    // - 这样可以在不改变 anchor 密度的前提下独立调节 sketch 稳定性。
    app.add_option("--sketch-kmer-size", opt.sketch_kmer_size,
                   "K-mer size used specifically for sketch construction.")
        ->default_val(21)
        ->check(CLI::Range(4, 31));
    // --batch-size：对齐阶段批大小。
    // 说明：
    // - 仅影响 alignSeq2Profile/alignSeq2Seq 的分批读取与并行粒度；
    // - 设为 0 时保持当前模式默认行为（不改变历史默认逻辑）。
    app.add_option("--batch-size", opt.batch_size,
                   "Alignment batch size (0 keeps mode default behavior).")
        ->default_val(0)
        ->check(CLI::Range(0, 100000000));

    app.add_option("-s,--score", opt.score_path,
                   "DNA5 alignment scoring matrix file for A/C/G/T/N (25 signed int8 scores).")
        ->check(CLI::ExistingFile);

    app.add_option("--gap-open", opt.gap_open,
                   "Gap open penalty used by reference alignment.")
        ->default_val(10)
        ->check(CLI::Range(0, 127));

    app.add_option("--gap-extend", opt.gap_extend,
                   "Gap extension penalty used by reference alignment.")
        ->default_val(2)
        ->check(CLI::Range(0, 127));

    // 开关参数：默认关闭，传入 --wfa 时设为 true
    app.add_flag("--wfa", opt.wfa,
        "Enable WFA alignment path (default: disabled).");

    // 比对模式开关：默认走 seq2profile，开启后切换为 seq2seq。
    app.add_flag("--seq2seq", opt.seq2seq,
        "Use seq2seq alignment pipeline instead of seq2profile (default: seq2profile).");


    app.add_flag("--keep-length", opt.keep_length,
        "Keep all reference sequences lengths unchanged. ");

    // workdir 管理：是否在完成后保留工作目录
    // --save-workdir：保留工作目录（默认会删除）。
    app.add_flag("--save-workdir", opt.save_workdir,
        "Keep the working directory after completion (default: remove). Useful for debugging.");
}

// logParsedOptions：把解析后的参数以漂亮的表格形式输出到日志
// 说明：此处的输出用于帮助用户和调试（打印被截断的长字符串、boolean 友好显示等），
// 不影响程序行为。若程序在无头环境运行（服务/容器），日志也便于审计和复现运行参数。
static void logParsedOptions(const Options& opt) {
    // Helper to convert values and truncate long strings for tidy display
    auto toString = [](const std::string& s, size_t maxLen) -> std::string {
        if (s.empty()) return "(empty)";
        if (s.size() <= maxLen) return s;
        return s.substr(0, maxLen - 3) + "...";
    };

    auto boolToStr = [](bool b) { return b ? "true" : "false"; };

    const size_t keyW = 14;
    const size_t valW = 60;
    const size_t innerW = keyW + 3 + valW; // "key : value"

    std::vector<std::pair<std::string, std::string>> rows = {
        {"input", toString(opt.input, valW)},
        {"output", toString(opt.output, valW)},
        {"workdir", toString(opt.workdir, valW)},
        {"ref", toString(opt.ref_path, valW)},
        {"ref-align", toString(opt.ref_align_path, valW)},
        {"msa_cmd", toString(opt.msa_cmd, valW)},
        {"threads", std::to_string(opt.threads)},
        {"kmer-size", std::to_string(opt.kmer_size)},
        {"sketch-kmer-size", std::to_string(opt.sketch_kmer_size)},
        {"kmer-window", std::to_string(opt.kmer_window)},
        {"cons_n", std::to_string(opt.cons_n)},
        {"sketch_size", std::to_string(opt.sketch_size)},
        {"batch-size", std::to_string(opt.batch_size)},
        {"score", toString(opt.score_path, valW)},
        {"gap-open", std::to_string(opt.gap_open)},
        {"gap-extend", std::to_string(opt.gap_extend)},
        {"wfa", boolToStr(opt.wfa)},
        {"seq2seq", boolToStr(opt.seq2seq)},
        {"keep-length", boolToStr(opt.keep_length)},
        {"save-workdir", boolToStr(opt.save_workdir)}
    };

    std::ostringstream oss;

    // top border
    oss << "+" << std::string(innerW, '-') << "+\n";

    // title centered
    const std::string title = " Parsed options ";
    size_t paddingLeft = 0;
    if (innerW > title.size()) paddingLeft = (innerW - title.size()) / 2;
    oss << "|" << std::string(paddingLeft, ' ') << title
        << std::string(innerW - paddingLeft - title.size(), ' ') << "|\n";

    // separator
    oss << "+" << std::string(innerW, '-') << "+\n";

    // rows
    for (auto &kv : rows) {
        oss << "| " << std::left << std::setw(keyW) << kv.first << " : "
            << std::setw(valW) << kv.second << "|\n";
    }

    // bottom border
    oss << "+" << std::string(innerW, '-') << "+";

    spdlog::info("\n{}", oss.str());
}

// ------------------------------------------------------------------
// CLI11 自定义格式器（美化选项输出）
// 说明：
// - 自定义 `make_option_opts` 可以在帮助中显示参数类型与默认值，便于用户理解；
// - 自定义 `make_usage` 提供更友好的使用示例与说明，方便新手快速上手。
// ------------------------------------------------------------------
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// 自定义参数展示样式（带默认值）
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";
		std::ostringstream out;
		out << " " << opt->get_type_name();
		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		return out.str();
	}

	// 提供使用示例；此处的 Example 可根据项目实际可执行名称更新
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./halign4 -i <ref.fa> -o <output.fa> -w </path/to/workdir> [options]\n\n"
			<< "Example:\n"
			<< "  ./halign4 -i ref.fa -o output.fa -w ./tmp -t 8\n\n";
		return out.str();
	}
};

// ------------------------------------------------------------------
// CLI11 自定义 validator：自动去除参数两侧空白
// 说明：有时用户在命令行中误加空格或复制粘贴带有换行，trim_whitespace 可以提高健壮性
// ------------------------------------------------------------------
inline CLI::Validator trim_whitespace = CLI::Validator(
	[](std::string& s) {
		auto start = s.find_first_not_of(" \t\n\r");
		auto end = s.find_last_not_of(" \t\n\r");
		s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
		return std::string();  // 空字符串表示验证通过
	}, ""
);



// ------------------------------------------------------------------
// 日志系统初始化
// 说明（中文注释）：
// - `setupLoggerWithFile(path)` 会创建一个异步的 spdlog 日志器，输出到控制台与指定目录下的日志文件；
// - 异步日志（async_logger）在高并发日志场景下能降低阻塞与 I/O 延迟，但需要在程序启动时配置线程池；
// - `setupLogger()` 仅输出到控制台，适合交互式调试或短期运行；
// - 两个函数都会设置默认日志级别（Debug/Info）并定期刷盘（flush_every），这有助于在崩溃时保留日志。
// ------------------------------------------------------------------
inline void setupLoggerWithFile(std::filesystem::path log_dir) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S] [%l] %v%$");

	std::filesystem::path log_file = log_dir / LOGGER_FILE;
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");

	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif
	spdlog::flush_every(std::chrono::seconds(3));
}

// 控制台日志（用于开发/调试）
inline void setupLogger() {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif

	spdlog::flush_every(std::chrono::seconds(3));
}

// 获取完整命令行字符串（用于日志记录和重现）
inline std::string getCommandLine(int argc, char** argv) {
	std::ostringstream cmd;
	for (int i = 0; i < argc; ++i) {
		cmd << argv[i];
		if (i != argc - 1) cmd << " ";
	}
	return cmd.str();
}

#endif // CONFIG_HPP

