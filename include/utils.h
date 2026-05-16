#ifndef HALIGN4_UTILS_H
#define HALIGN4_UTILS_H

// ================================================================
// utils.h - 文件与序列 I/O 公共声明（详细中文注释）
//
// 本头文件包含两个主要功能块：
//  1) file_io: 一组与文件/路径/网络获取相关的工具函数（copy/download/fetch/prep目录等），
//     目的在于把所有与文件系统交互的逻辑集中管理，便于统一错误处理与测试。
//  2) seq_io: 一组基于 kseq 的 FASTA/FASTQ 读写抽象（KseqReader / FastaWriter / SeqRecord），
//     并提供序列清洗函数（cleanSequence）和便捷的 reader/writer 构造器。
// 此处只声明接口；具体实现位于 src/utils/file_io.cpp 和 src/utils/seq_io.cpp 中。
// ================================================================

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <regex>
#include <array>
#include <cstdint>
#include <unordered_map>

#if __has_include(<curl/curl.h>)
#include <curl/curl.h>
#define HALIGN4_HAVE_LIBCURL 1
#endif

namespace fs = std::filesystem;
using FilePath = std::filesystem::path;


namespace file_io
{
    // 下面是 file_io 中提供的函数声明，并在每个函数前添加中文说明，帮助使用者理解语义与边界。

    // 格式化文件系统错误信息为可读字符串（带路径与错误消息）
    // 用法：在抛异常或记录日志时调用，统一错误输出格式
    std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec);

    // 要求路径存在（任何类型）；若不存在或检查失败则抛出 runtime_error
    void requireExists(const FilePath& p, std::string_view what);

    // 要求路径是常规文件；若不是则抛出 runtime_error
    void requireRegularFile(const FilePath& p, std::string_view what);

    // 要求路径是目录；若不是则抛出 runtime_error
    void requireDirectory(const FilePath& p, std::string_view what);

    // 确保目录存在；若不存在则创建（包括父目录）。出错抛异常。
    void ensureDirectoryExists(const FilePath& p, std::string_view what = "directory");

    // 判断路径是否为空（对文件适用）
    bool isEmpty(const FilePath& p);

    // 核心：准备工作目录（不存在则创建；存在则要求为目录；可选要求为空）
    // - workdir: 需准备的目录路径
    // - must_be_empty: 若 true，则目录存在时必须为空，否则抛错
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty);

    // 可选：确保输出文件的父目录存在（写文件前调用）
    void ensureParentDirExists(const FilePath& out_file);

    // 判断给定的 FilePath 是否表示一个远程 URL（例如 http://、https://、ftp://、s3:// 等）
    // 返回 true 表示应当通过网络下载（远程链接）；返回 false 表示为本地路径。
    // 注意：该判断为启发式；对 file:// 或其他边缘 scheme，调用方应根据需要额外处理。
    bool isUrl(const FilePath& p);

    // copyFile：复制文件并覆盖目标；优先使用 std::filesystem::copy_file。
    // 若遇到跨设备错误（EXDEV），回退到基于流的拷贝以确保可跨分区复制。
    // 出错抛异常。
    void copyFile(const FilePath& src, const FilePath& dst);

    // downloadFile：把远程 URL 下载到本地路径 dst
    // - 如果项目编译时检测到 libcurl，则使用其 API（更可靠、可控）
    // - 否则回退到命令行工具（curl/wget），此方式依赖运行时环境
    // 出错抛异常并尝试删除局部残留文件。
    void downloadFile(const std::string& url, const FilePath& dst);

    // fetchFile：如果 srcOrUrl 是 URL，则调用 downloadFile，否则调用 copyFile
    // 这是上层调用统一接口，屏蔽本地/远程差异
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst);

    // 递归删除路径（文件或目录），失败时抛出异常。慎用！
    void removeAll(const FilePath& p);

    // 读取整个小文件到 std::string（仅适用于文件体积适中时）
    std::string readFileToString(const FilePath& p);

}


namespace seq_io
{
#if defined(FilePath)
    using FilePath = ::FilePath;
#else
    using FilePath = std::filesystem::path;
#endif

    // SeqRecord：用于在内存中表示一条序列记录（FASTA/FASTQ）
    // 字段语义：
    //  - id: header 名称，通常是 '>' 后第一个空格前的 token
    //  - desc: header 剩余描述信息（可选）
    //  - seq: 序列字符串（A/C/G/T/U/-/N 等）
    //  - qual: FASTQ 质量字符串（若为 FASTQ，可填充；FASTA 情况可留空）
    struct SeqRecord
    {
        std::string id;     // header name（到空格前）
        std::string desc;   // header 其余部分
        std::string seq;    // 序列
        std::string qual;   // FASTQ 质量字符串（可选），实现中需要访问该字段
        int n_num = 0;      // 序列中 N/n 的数量
    };
    using SeqRecords     = std::vector<seq_io::SeqRecord>;
	using SeqRecordMap   = std::unordered_map<std::string, SeqRecord>;

    // ------------------------------------------------------------------
    // SamRecord：用于在内存中表示一条 SAM 比对记录
    // 字段语义（遵循 SAM 格式标准）：
    //  - qname: Query 名称（QNAME 列）
    //  - flag: SAM flag 位标志（FLAG 列，例如 0=正向未配对，16=反向链等）
    //  - rname: Reference 序列名称（RNAME 列，"*"表示未比对）
    //  - pos: 1-based 比对起始位置（POS 列，0 或 "*" 表示未知）
    //  - mapq: Mapping quality（MAPQ 列，0-255，255 表示未知）
    //  - cigar: CIGAR 字符串（CIGAR 列，例如 "100M5I95M"，"*"表示未知）
    //  - rnext: Mate/next reference 名称（RNEXT 列，"="表示相同，"*"表示未知）
    //  - pnext: Mate 1-based 位置（PNEXT 列）
    //  - tlen: Template length（TLEN 列，有符号整数，0 表示未知）
    //  - seq: Query 序列（SEQ 列）
    //  - qual: Query 质量字符串（QUAL 列，"*"表示未知）
    //  - opt: 可选 TAG 字段（OPT 列，从第 12 列开始的所有内容）
    //
    // 设计说明：
    // - 所有字符串字段使用 std::string（拥有字段所有权，无生命周期问题）
    // - 适用于 SamReader/SamWriter 的读写和长期存储
    // - 字段默认值遵循 SAM 规范："*" 表示未知/缺失
    // - 性能说明：使用 std::string 会有字符串拷贝开销，但换来更简单的生命周期管理
    // ------------------------------------------------------------------
    struct SamRecord
    {
        std::string qname;
        std::uint16_t flag = 0;
        std::string rname = "*";
        std::uint32_t pos = 0;     // 1-based; 0 means '*'
        std::uint8_t mapq = 0;
        std::string cigar = "*";
        std::string rnext = "*";
        std::uint32_t pnext = 0;
        std::int32_t tlen = 0;
        std::string seq = "*";
        std::string qual = "*";
        std::string opt;  // 可选 TAG 字段（包含前导 '\t' 或不包含都可）
    };

    // makeCleanTable / clean_table:
    // - 用于把任意字符映射到一个受限的碱基字符集合（'A','C','G','T','U','N','-'），
    // - 实现为 256 大小的查表（针对 unsigned char），避免运行时的分支/条件判断，极大提升批量清洗性能；
    // - cleanSequence 提供了 in-place 清洗接口，对大序列重复使用时能减少临时分配。
    constexpr std::array<std::uint8_t, 256> makeCleanTable()
    {
        std::array<std::uint8_t, 256> table{};

        for (std::size_t i = 0; i < table.size(); ++i) {
            table[i] = static_cast<std::uint8_t>('N');
        }

        table[static_cast<unsigned char>('A')] = 'A';
        table[static_cast<unsigned char>('a')] = 'A';
        table[static_cast<unsigned char>('C')] = 'C';
        table[static_cast<unsigned char>('c')] = 'C';
        table[static_cast<unsigned char>('G')] = 'G';
        table[static_cast<unsigned char>('g')] = 'G';
        table[static_cast<unsigned char>('T')] = 'T';
        table[static_cast<unsigned char>('t')] = 'T';
        table[static_cast<unsigned char>('U')] = 'U';
        table[static_cast<unsigned char>('u')] = 'U';
        table[static_cast<unsigned char>('-')] = '-';
    	table[static_cast<unsigned char>('.')] = '-';

        return table;
    }

    // 头文件内的“单一定义”常量表：inline constexpr 是推荐写法
    inline constexpr auto clean_table = makeCleanTable();
    inline void cleanSequence(std::string& seq)
    {
        for (char& ch : seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
        }
    }
    inline void cleanSequence(SeqRecord& seq)
    {
        seq.n_num = 0;
        for (char& ch : seq.seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
            if (ch == 'N') {
                ++seq.n_num;
            }
        }
    }

    // ISequenceReader：抽象读取器接口，支持 next(SeqRecord&) 返回 false 表示 EOF
    class ISequenceReader
    {
    public:
        virtual ~ISequenceReader() = default;
        virtual bool next(SeqRecord& rec) = 0; // false => EOF
    };

    // KseqReader：基于 kseq 实现的高效序列读取器
    // 说明：实现位于 src/utils/seq_io.cpp，构造函数会根据是否启用 zlib 选择 gzopen/gzread 或 fopen/fread
    class KseqReader final : public ISequenceReader
    {
    public:
        explicit KseqReader(const FilePath& file_path);
        ~KseqReader() override;

        KseqReader(const KseqReader&) = delete;
        KseqReader& operator=(const KseqReader&) = delete;

        KseqReader(KseqReader&&) noexcept;
        KseqReader& operator=(KseqReader&&) noexcept;

        bool next(SeqRecord& rec) override;

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // SeqWriter：通用序列/比对结果写出器
    //
    // 支持格式：
    // - FASTA：写出 SeqRecord（id/desc/seq）
    // - SAM  ：写出最小可用的 SAM（header + alignment records）
    //
    // 说明：
    // - 默认使用 8MiB 的内部聚合缓冲，尽可能减少磁盘写调用次数；
    // - 为了兼容旧用法，保留 write(const SeqRecord&) 作为 FASTA 写出的默认入口。
    class SeqWriter
    {
    public:
        enum class Format : std::uint8_t
        {
            fasta = 0,
            sam   = 1,
        };

        // FASTA（默认）
        explicit SeqWriter(const FilePath& file_path, std::size_t line_width = 80);
        explicit SeqWriter(const FilePath& file_path, Format fmt, std::size_t line_width, std::size_t buffer_threshold_bytes);

        SeqWriter(const FilePath& file_path, std::size_t line_width, std::size_t buffer_threshold_bytes);

        SeqWriter(const SeqWriter&) = delete;
        SeqWriter& operator=(const SeqWriter&) = delete;

        SeqWriter(SeqWriter&&) noexcept = default;
        SeqWriter& operator=(SeqWriter&&) noexcept = default;


        // SAM 构造（工厂函数）
        static SeqWriter Sam(const FilePath& file_path, std::size_t buffer_threshold_bytes = 8ULL * 1024ULL * 1024ULL);

        // 当前模式
        Format format() const noexcept { return format_; }

        // ------------------------- FASTA -------------------------
        void writeFasta(const SeqRecord& rec);

        // 兼容：默认 write() 写 FASTA
        void write(const SeqRecord& rec) { writeFasta(rec); }

        // ------------------------- SAM -------------------------
        void writeSamHeader(std::string_view header_text);

        void writeSam(const SamRecord& r);

        // flush():
        // - 先 flush 内部聚合缓冲，再 flush ofstream
        void flush();

        ~SeqWriter();

    private:

        std::ofstream out_;
        Format format_{Format::fasta};
        std::size_t line_width_{80};

        std::string buffer_;
        std::size_t buffer_threshold_bytes_{8ULL * 1024ULL * 1024ULL};

        bool sam_header_written_{false};

        void flushBuffer_();
        void appendOrFlush_(std::string_view s);
        static void appendWrapped_(std::string& dst, std::string_view s, std::size_t width);
    };

    inline std::unique_ptr<ISequenceReader> openKseqReader(const FilePath& file_path)
    {
        return std::make_unique<KseqReader>(file_path);
    }

    // ------------------------------------------------------------------
    // 类：SamReader
    // 功能：读取 SAM 文件并解析为 SamRecord
    //
    // 说明：
    // 1. SAM 格式为制表符分隔的文本文件，每行包含一条比对记录
    // 2. SAM 文件可包含以 '@' 开头的 header 行，本类会自动跳过
    // 3. 读取时解析所有 SAM 必需字段（QNAME、FLAG、RNAME、POS、MAPQ、CIGAR、SEQ、QUAL 等）
    // 4. 支持大文件高效读取（使用大缓冲区优化）
    //
    // 性能说明：
    // - 使用 std::getline 逐行读取，配合大缓冲区减少系统调用
    // - 使用 string_view 进行字段分割，避免不必要的字符串拷贝
    // - 默认缓冲区大小 8MiB，适合大规模 SAM 文件
    // - 内部保存字段字符串以确保 string_view 生命周期安全
    // ------------------------------------------------------------------
    class SamReader
    {
    public:
        // 构造函数：打开 SAM 文件并初始化读取器
        // 参数：
        //   - file_path: SAM 文件路径
        //   - buffer_size: 输入缓冲区大小（默认 8MiB）
        explicit SamReader(const FilePath& file_path, std::size_t buffer_size = 8ULL * 1024ULL * 1024ULL);

        ~SamReader();

        // 禁止拷贝，允许移动
        SamReader(const SamReader&) = delete;
        SamReader& operator=(const SamReader&) = delete;
        SamReader(SamReader&&) noexcept;
        SamReader& operator=(SamReader&&) noexcept;


        bool next(SamRecord& rec);

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // ------------------------------------------------------------------
    // 函数：convertSamToFasta
    // 功能：读取 SAM 文件并转换为 FASTA 文件
    //
    // 说明：
    // 1. 从 SAM 文件中提取 QNAME 和 SEQ 字段，写入 FASTA 格式
    // 2. 自动跳过 SAM header 行（以 '@' 开头）
    // 3. 保留 query 序列的名称和序列内容，丢弃比对信息和质量值
    // 4. 支持大文件高效转换（使用缓冲区优化）
    //
    // 参数：
    //   - sam_path: 输入 SAM 文件路径
    //   - fasta_path: 输出 FASTA 文件路径
    //   - line_width: FASTA 每行宽度（默认 80，0 表示不换行）
    //
    // 性能：
    // - 使用 8MiB 输入/输出缓冲区，减少 I/O 系统调用
    // - 流式处理，内存占用与文件大小无关
    //
    // 异常：文件打开失败或解析错误时抛出 std::runtime_error
    // ------------------------------------------------------------------
    void convertSamToFasta(const FilePath& sam_path, const FilePath& fasta_path, std::size_t line_width = 80);


    // ------------------------------------------------------------------
    // 函数：makeSamRecord
    // 功能：从 SeqRecord 和比对信息构建 SAM 记录
    //
    // 参数：
    //   - query: 查询序列的 SeqRecord
    //   - ref_name: 参考序列名称
    //   - cigar_str: CIGAR 字符串（例如 "100M5I95M"）
    //   - pos: 1-based 比对起始位置（默认 1）
    //   - mapq: mapping quality（默认 60）
    //   - flag: SAM flag（默认 0）
    //
    // 返回：SamRecord（seq_io 命名空间下的独立类型）
    //
    // 说明：
    //   1. 这是一个便捷函数，封装了从 SeqRecord 到 SamRecord 的转换逻辑
    //   2. 质量值处理：若 query.qual 为空，则使用 SAM 默认值 "*"
    //   3. 所有 string_view 字段的生命周期由调用者保证（query 和 cigar_str 必须在使用期间有效）
    //   4. 简化的 SAM 记录，后续可扩展添加更多字段（AS、NM 等）
    //
    // 性能：O(1)，仅赋值操作，无内存分配
    // ------------------------------------------------------------------
    SamRecord makeSamRecord(
        const SeqRecord& query,
        std::string_view ref_name,
        std::string_view cigar_str,
        std::uint32_t pos = 1,
        std::uint8_t mapq = 60,
        std::uint16_t flag = 0
    );

    // ------------------------------------------------------------------
    // 函数：samRecordToSeqRecord
    // 功能：把 SamRecord 中的序列信息（QNAME/SEQ/QUAL）转换为 SeqRecord
    //
    // 为什么需要它：
    // - SamReader::next() 返回的是 SamRecord（字段多，包含所有 SAM 列）
    // - 但很多后续处理（例如写 FASTA、cleanSequence、与原有代码对接）更习惯使用 SeqRecord
    // - 这个函数把"从 SamRecord 提取序列"的重复逻辑集中到一处，避免到处写拷贝代码
    //
    // 输入：
    // - sam_rec: SamRecord（通常来自 SamReader::next），其中 qname/seq/qual 为 std::string
    // - keep_qual: 是否把 qual 拷贝到 SeqRecord（FASTA 场景通常不需要 qual，建议 false）
    //
    // 输出：
    // - 返回一个 SeqRecord：
    //   - id   = sam_rec.qname
    //   - desc = ""（SAM 没有 FASTA header 的 desc 概念）
    //   - seq  = sam_rec.seq
    //   - qual = (keep_qual ? sam_rec.qual : "")
    //
    // 性能与正确性：
    // - 由于 SamRecord 已经拥有字符串的所有权，这里会发生字符串拷贝
    // - 复杂度 O(|qname| + |seq| + |qual|)，对大文件通常 IO 才是瓶颈
    // ------------------------------------------------------------------
    inline SeqRecord samRecordToSeqRecord(const SamRecord& sam_rec, bool keep_qual = false)
    {
        SeqRecord rec;
        rec.id = sam_rec.qname;
        rec.desc.clear();
        rec.seq = sam_rec.seq;

        if (keep_qual) {
            // SAM 的 qual 可能是 "*"，调用方如需可自行再做过滤。
            rec.qual = sam_rec.qual;
        } else {
            rec.qual.clear();
        }

        // 统计当前序列中 N/n 的数量，保证不同构造路径字段语义一致。
        rec.n_num = 0;

        return rec;
    }

} // namespace seq_io

namespace cmd
{
    struct BuildOptions
    {
        bool quiet = true;                  // 是否追加静默重定向
        bool close_stdin = true;            // 是否追加 "< /dev/null"
        bool detect_stdout_redirect = true; // 是否根据 '>' 决定静默策略
    };

    // cmd_template 必须包含 {input} 和 {output}；可选包含 {thread}
    // Linux-only 规则：
    // - 若模板中检测到 '>'（stdout 已被重定向）：追加 "2>/dev/null"
    // - 否则追加 "> /dev/null 2>&1"
    // - 可选追加 "< /dev/null" 关闭 stdin
    std::string buildCommand(std::string cmd_template,
                             const std::string& input_path,
                             const std::string& output_path,
                             int thread,
                             const BuildOptions& opt = BuildOptions{});

    // 执行命令并返回真实退出码：
    // - 正常退出：返回 exit code
    // - 信号终止：返回 128 + signal
    // - system() 自身失败：返回 -1
    int runCommand(const std::string& command);

    // 用 tiny.fasta 做自检，运行模板命令并检查输出文件是否生成
    bool testCommandTemplate(const std::string& cmd_template, const FilePath& workdir, int thread = 1);

} // namespace cmd


// ================================================================
// ProgressBar - 进度条工具类（单行刷新，无需预知总数）
//
// 设计目标：
// 1) 统一项目中所有进度条的实现，避免代码冗余
// 2) 支持不预先知道总数的场景（只显示已处理数 + 速率）
// 3) 线程安全考虑：本类本身不加锁，适用于串行区调用
//    若需在并行区使用，应由调用方保证串行访问（如 OpenMP 串行区）
//
// 使用示例：
//   ProgressBar progress("align");
//   for (...) {
//       // ... 处理工作 ...
//       progress.tick();  // 处理一条
//   }
//   progress.done();  // 强制输出并换行
// ================================================================
class ProgressBar {
public:
    // 构造函数
    // - label: 进度条前缀标签（如 "preprocess"、"align"、"merge"）
    // - report_interval: 每处理多少条刷新一次（默认 1000）
    explicit ProgressBar(std::string label, std::size_t report_interval = 1000) noexcept
        : label_(std::move(label))
        , report_interval_(report_interval)
        , count_(0)
        , next_report_(report_interval)
        , start_time_(std::chrono::steady_clock::now())
    {}

    // 禁用拷贝（避免意外复制导致的状态混乱）
    ProgressBar(const ProgressBar&) = delete;
    ProgressBar& operator=(const ProgressBar&) = delete;

    // 允许移动
    ProgressBar(ProgressBar&&) noexcept = default;
    ProgressBar& operator=(ProgressBar&&) noexcept = default;

    // tick：处理完一条记录后调用，内部自动判断是否需要刷新
    void tick() noexcept {
        ++count_;
        if (count_ >= next_report_) {
            print(false);
            next_report_ = ((count_ / report_interval_) + 1) * report_interval_;
        }
    }

    // tick：处理完多条记录后调用（批量更新，减少调用开销）
    void tick(std::size_t n) noexcept {
        count_ += n;
        if (count_ >= next_report_) {
            print(false);
            next_report_ = ((count_ / report_interval_) + 1) * report_interval_;
        }
    }

    // done：强制输出最终状态并换行（通常在循环结束后调用）
    void done() noexcept {
        print(true);
        std::fprintf(stderr, "\n");
    }

    // 获取当前已处理的数量
    [[nodiscard]] std::size_t count() const noexcept { return count_; }

    // 重置进度条（用于复用同一实例处理多个阶段）
    void reset() noexcept {
        count_ = 0;
        next_report_ = report_interval_;
        start_time_ = std::chrono::steady_clock::now();
    }

    // 重置进度条并更换标签
    void reset(std::string new_label) noexcept {
        label_ = std::move(new_label);
        reset();
    }

private:
    // print：输出进度信息到 stderr（使用 \r 回到行首覆盖）
    // - force: 若为 true 则强制输出，否则仅在达到阈值时输出
    void print(bool force) const noexcept {
        if (!force && count_ < next_report_) {
            return;
        }

        const auto now = std::chrono::steady_clock::now();
        const double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(
            now - start_time_).count();
        const double rate = (elapsed_sec > 1e-6)
            ? (static_cast<double>(count_) / elapsed_sec)
            : 0.0;

        // ANSI 颜色码：绿色 = \033[32m，重置 = \033[0m
        // 末尾空格用于覆盖可能残留的旧字符
        std::fprintf(stderr,
                     "\r\033[32m[%s] processed=%zu  elapsed=%.1fs  rate=%.0f seq/s\033[0m   ",
                     label_.c_str(), count_, elapsed_sec, rate);
        std::fflush(stderr);
    }

    std::string label_;                                  // 进度条前缀标签
    std::size_t report_interval_;                        // 刷新间隔（每处理多少条刷新一次）
    std::size_t count_;                                  // 已处理记录数
    std::size_t next_report_;                            // 下次刷新阈值
    std::chrono::steady_clock::time_point start_time_;   // 开始时间
};

#endif //HALIGN4_UTILS_H

