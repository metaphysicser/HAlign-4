// 序列文件 I/O 模块：FASTA/FASTQ/SAM 读写
// - KseqReader：高性能 FASTA/FASTQ 读取器
// - SeqWriter：FASTA/SAM 写入器（支持批量缓冲）
// - SamReader：SAM 格式读取器
// - 辅助函数：格式转换、记录构建

#include "utils.h"
#include <cstdio>   // fopen, fclose, fread, setvbuf
#include <cstdlib>  // malloc, free

#ifdef _DEBUG
#include <spdlog/spdlog.h>
#endif

// zlib 支持：可选的 gzip 压缩文件读取
#if __has_include(<zlib.h>)
    #include <zlib.h>
    #define HALIGN4_HAVE_ZLIB 1
#else
    #define HALIGN4_HAVE_ZLIB 0
#endif

#include "kseq.h"

// kseq 初始化：根据是否有 zlib 选择不同的 read 函数
#if HALIGN4_HAVE_ZLIB
    KSEQ_INIT(gzFile, gzread)
#else
    static int fileRead(std::FILE* fp, void* buf, int len)
    {
        const std::size_t n = std::fread(buf, 1, static_cast<std::size_t>(len), fp);
        return static_cast<int>(n);
    }
    KSEQ_INIT(std::FILE*, fileRead)
#endif

namespace seq_io
{
    // KseqReader::Impl 结构体：封装底层文件描述符和 kseq 对象
    struct KseqReader::Impl
    {
#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};
#else
        std::FILE* fp{nullptr};
#endif
        kseq_t* seq{nullptr};
        FilePath file_path;
        char* io_buf{nullptr};
        std::size_t io_buf_size{8 << 20};  // 默认 8 MiB
    };

    // 错误处理：构建包含文件路径的异常对象
    static std::runtime_error makeIoError(const std::string& msg, const FilePath& p)
    {
        return std::runtime_error(msg + ": " + p.string());
    }

    // 字符串转换：将 kstring_t 转换为 std::string
    static void assignKstring(std::string& dst, const kstring_t& ks)
    {
        if (ks.s && ks.l > 0) dst.assign(ks.s, ks.l);
        else dst.clear();
    }

    // 构建完整 FASTA/FASTQ header。
    // kseq 会把 header 的第一个空白前内容放在 name，后面的注释放在 comment。
    // 为避免只保存 name 造成大量重复 header，这里把 name + comment 合并保存。
    static void assignKseqFullHeader(std::string& dst, const kstring_t& name, const kstring_t& comment)
    {
        dst.clear();

        if (name.s && name.l > 0) {
            dst.assign(name.s, name.l);
        }

        if (comment.s && comment.l > 0) {
            if (!dst.empty()) {
                dst.push_back(' ');
            }
            dst.append(comment.s, comment.l);
        }
    }

    // KseqReader 构造函数
    KseqReader::KseqReader(const FilePath& file_path)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path = file_path;

#if HALIGN4_HAVE_ZLIB
        impl_->fp = gzopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }
        gzbuffer(impl_->fp, static_cast<int>(impl_->io_buf_size));
        impl_->seq = kseq_init(impl_->fp);
#else
        impl_->fp = std::fopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 设置大缓冲区以加速 I/O
        impl_->io_buf = static_cast<char*>(std::malloc(impl_->io_buf_size));
        if (impl_->io_buf) {
            if (setvbuf(impl_->fp, impl_->io_buf, _IOFBF, static_cast<size_t>(impl_->io_buf_size)) != 0) {
                std::free(impl_->io_buf);
                impl_->io_buf = nullptr;
            }
        }

        impl_->seq = kseq_init(impl_->fp);
#endif

        if (!impl_->seq) {
            throw makeIoError("failed to init kseq", file_path);
        }
    }

    // KseqReader 析构函数：按正确顺序释放资源
    KseqReader::~KseqReader()
    {
        if (!impl_) return;

        if (impl_->seq) {
            kseq_destroy(impl_->seq);
            impl_->seq = nullptr;
        }

#if HALIGN4_HAVE_ZLIB
        if (impl_->fp) {
            gzclose(impl_->fp);
            impl_->fp = nullptr;
        }
#else
        if (impl_->fp) {
            std::fclose(impl_->fp);
            impl_->fp = nullptr;
        }
        if (impl_->io_buf) {
            std::free(impl_->io_buf);
            impl_->io_buf = nullptr;
        }
#endif
    }

    // 移动语义
    KseqReader::KseqReader(KseqReader&& other) noexcept = default;
    KseqReader& KseqReader::operator=(KseqReader&& other) noexcept = default;

    // 读取下一条 FASTA/FASTQ 记录
    bool KseqReader::next(SeqRecord& rec)
    {
        if (!impl_ || !impl_->seq) {
            throw std::runtime_error("KseqReader is not initialized");
        }

        const int ret = kseq_read(impl_->seq);

        if (ret >= 0) {
            // 保存完整 header，而不是只保存 kseq 的 name 字段。
            // 原始 header: >name comment
            // kseq 拆分为: name=name, comment=comment
            // 这里将二者合并到 rec.id 中，保证后续只使用 rec.id 的流程也不会丢失注释部分。
            assignKseqFullHeader(rec.id, impl_->seq->name, impl_->seq->comment);

            // comment 已经合并进 rec.id。清空 rec.desc，避免 SeqWriter::writeFasta 再次追加，
            // 导致输出 header 变成 ">name comment comment"。
            rec.desc.clear();

            assignKstring(rec.seq, impl_->seq->seq);
            rec.n_num = 0;
            return true;
        }

        if (ret == -1) {
            return false;  // EOF
        }

        throw std::runtime_error("kseq_read() failed with code " + std::to_string(ret) +
                                 " for file: " + impl_->file_path.string());
    }

    // SeqWriter 构造函数
    SeqWriter::SeqWriter(const FilePath& file_path, std::size_t line_width)
        : SeqWriter(file_path, Format::fasta, line_width, 8ULL * 1024ULL * 1024ULL)
    {}

    SeqWriter::SeqWriter(const FilePath& file_path, std::size_t line_width, std::size_t buffer_threshold_bytes)
        : SeqWriter(file_path, Format::fasta, line_width, buffer_threshold_bytes)
    {}

    SeqWriter::SeqWriter(const FilePath& file_path, Format fmt, std::size_t line_width, std::size_t buffer_threshold_bytes)
        : out_(file_path, std::ios::binary),
          format_(fmt),
          line_width_(line_width == 0 ? 80 : line_width),
          buffer_threshold_bytes_(buffer_threshold_bytes)
    {
        if (!out_) {
            throw makeIoError("failed to open output", file_path);
        }

        if (buffer_threshold_bytes_ > 0) {
            buffer_.reserve(std::min<std::size_t>(buffer_threshold_bytes_, 8ULL * 1024ULL * 1024ULL));
        }
    }

    SeqWriter SeqWriter::Sam(const FilePath& file_path, std::size_t buffer_threshold_bytes)
    {
        return SeqWriter(file_path, Format::sam, /*line_width=*/0, buffer_threshold_bytes);
    }

    SeqWriter::~SeqWriter()
    {
        try {
            flush();
        } catch (...) {
            // best-effort
        }
    }

    // 刷新缓冲区并释放容量
    void SeqWriter::flushBuffer_()
    {
        if (buffer_.empty()) return;
        out_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        buffer_.clear();
        buffer_.shrink_to_fit();  // 释放未使用容量
    }

    // 追加内容到缓冲区或直接写入
    void SeqWriter::appendOrFlush_(std::string_view s)
    {
        if (!out_) throw std::runtime_error("SeqWriter output stream is not ready");

        if (buffer_threshold_bytes_ == 0) {
            out_.write(s.data(), static_cast<std::streamsize>(s.size()));
            return;
        }

        buffer_.append(s.data(), s.size());
        if (buffer_.size() >= buffer_threshold_bytes_) {
            flushBuffer_();
        }
    }

    // 追加内容并支持行折叠
    void SeqWriter::appendWrapped_(std::string& dst, std::string_view s, std::size_t width)
    {
        if (width == 0) {
            dst.append(s.data(), s.size());
            return;
        }

        for (std::size_t i = 0; i < s.size(); i += width) {
            const std::size_t n = (i + width <= s.size()) ? width : (s.size() - i);
            dst.append(s.data() + i, n);
            dst.push_back('\n');
        }
    }

    // 写入 FASTA 记录
    void SeqWriter::writeFasta(const SeqRecord& rec)
    {
        if (format_ != Format::fasta) {
            throw std::runtime_error("SeqWriter::writeFasta called but writer is not in FASTA mode");
        }

        const std::size_t width = (line_width_ == 0) ? 80 : line_width_;

        std::string recordbuf;
        recordbuf.reserve(1 + rec.id.size() + (rec.desc.empty() ? 0 : 1 + rec.desc.size()) + 1 +
                          rec.seq.size() + (rec.seq.size() / width) + 2);

        // header
        recordbuf.push_back('>');
        recordbuf.append(rec.id);
        if (!rec.desc.empty()) {
            recordbuf.push_back(' ');
            recordbuf.append(rec.desc);
        }
        recordbuf.push_back('\n');

        // sequence (wrapped)
        if (rec.seq.empty()) {
            recordbuf.push_back('\n');
        } else {
            for (std::size_t i = 0; i < rec.seq.size(); i += width) {
                const std::size_t n = (i + width <= rec.seq.size()) ? width : (rec.seq.size() - i);
                recordbuf.append(rec.seq.data() + i, n);
                recordbuf.push_back('\n');
            }
        }

        appendOrFlush_(recordbuf);
    }

    // 写入 SAM header
    void SeqWriter::writeSamHeader(std::string_view header_text)
    {
        if (format_ != Format::sam) {
            throw std::runtime_error("SeqWriter::writeSamHeader called but writer is not in SAM mode");
        }

        if (header_text.empty()) return;

        appendOrFlush_(header_text);
        if (header_text.back() != '\n') {
            appendOrFlush_("\n");
        }

        sam_header_written_ = true;
    }

    // 写入 SAM 记录
    void SeqWriter::writeSam(const SamRecord& r)
    {
        if (format_ != Format::sam) {
            throw std::runtime_error("SeqWriter::writeSam called but writer is not in SAM mode");
        }

        std::string line;
        line.reserve(r.qname.size() + r.rname.size() + r.cigar.size() + r.rnext.size() + r.seq.size() + r.qual.size() + 64);

        line.append(r.qname.data(), r.qname.size());
        line.push_back('\t');
        line.append(std::to_string(r.flag));
        line.push_back('\t');
        line.append(r.rname.data(), r.rname.size());
        line.push_back('\t');
        line.append(std::to_string(r.pos));
        line.push_back('\t');
        line.append(std::to_string(static_cast<unsigned int>(r.mapq)));
        line.push_back('\t');
        line.append(r.cigar.data(), r.cigar.size());
        line.push_back('\t');
        line.append(r.rnext.data(), r.rnext.size());
        line.push_back('\t');
        line.append(std::to_string(r.pnext));
        line.push_back('\t');
        line.append(std::to_string(r.tlen));
        line.push_back('\t');
        line.append(r.seq.data(), r.seq.size());
        line.push_back('\t');
        line.append(r.qual.data(), r.qual.size());

        if (!r.opt.empty()) {
            if (r.opt.front() != '\t') line.push_back('\t');
            line.append(r.opt.data(), r.opt.size());
        }

        line.push_back('\n');
        appendOrFlush_(line);
    }

    // 刷新输出
    void SeqWriter::flush()
    {
        if (!out_) return;
        flushBuffer_();
        out_.flush();
    }

    // 构建 SAM 记录
    SamRecord makeSamRecord(
        const SeqRecord& query,
        std::string_view ref_name,
        std::string_view cigar_str,
        std::uint32_t pos,
        std::uint8_t mapq,
        std::uint16_t flag)
    {
        SamRecord sam_rec;
        sam_rec.qname = query.id;
        sam_rec.flag  = flag;
        sam_rec.rname.assign(ref_name.data(), ref_name.size());
        sam_rec.pos   = pos;
        sam_rec.mapq  = mapq;
        sam_rec.cigar.assign(cigar_str.data(), cigar_str.size());
        sam_rec.seq   = query.seq;

        if (!query.qual.empty()) {
            sam_rec.qual = query.qual;
        }

        return sam_rec;
    }

    // SamReader::Impl 结构体
    struct SamReader::Impl
    {
        std::ifstream in_;
        FilePath file_path_;
        std::string line_buffer_;
        char* io_buf_{nullptr};
        std::size_t io_buf_size_{0};

        // 存储字段副本，避免 string_view 悬空
        std::string qname_storage_;
        std::string rname_storage_;
        std::string cigar_storage_;
        std::string rnext_storage_;
        std::string seq_storage_;
        std::string qual_storage_;
        std::string opt_storage_;
    };

    // SamReader 构造函数
    SamReader::SamReader(const FilePath& file_path, std::size_t buffer_size)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path_ = file_path;
        impl_->io_buf_size_ = buffer_size;

        impl_->in_.open(file_path, std::ios::in);
        if (!impl_->in_) {
            throw makeIoError("failed to open SAM file", file_path);
        }

        // 提升输入缓冲区大小
        if (buffer_size > 0) {
            impl_->io_buf_ = static_cast<char*>(std::malloc(buffer_size));
            if (impl_->io_buf_) {
                impl_->in_.rdbuf()->pubsetbuf(impl_->io_buf_, static_cast<std::streamsize>(buffer_size));
            }
        }

        impl_->line_buffer_.reserve(4096);
    }

    // SamReader 析构函数
    SamReader::~SamReader()
    {
        if (!impl_) return;

        if (impl_->in_.is_open()) {
            impl_->in_.close();
        }

        if (impl_->io_buf_) {
            std::free(impl_->io_buf_);
            impl_->io_buf_ = nullptr;
        }
    }

    // 移动语义
    SamReader::SamReader(SamReader&& other) noexcept = default;
    SamReader& SamReader::operator=(SamReader&& other) noexcept = default;

    // 读取下一条 SAM 记录
    bool SamReader::next(SamRecord& rec)
    {
        if (!impl_) {
            throw std::runtime_error("SamReader is not initialized");
        }

        while (std::getline(impl_->in_, impl_->line_buffer_)) {
            const std::string_view line(impl_->line_buffer_);

            // 跳过空行和 header 行
            if (line.empty() || line[0] == '@') {
                continue;
            }

            // 解析 SAM 必需字段 (11 列)
            std::array<std::string_view, 11> f{};
            std::size_t field_idx = 0;
            std::size_t start = 0;

            for (std::size_t i = 0; i <= line.size() && field_idx < f.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    f[field_idx] = line.substr(start, i - start);
                    start = i + 1;
                    ++field_idx;
                }
            }

            if (field_idx < f.size()) {
                throw std::runtime_error(
                    "invalid SAM record (missing required fields): " + impl_->file_path_.string());
            }

            // 可选字段
            std::string_view opt_fields;
            if (start < line.size()) {
                opt_fields = line.substr(start);
            }

            // 填充 SamRecord
            rec.qname.assign(f[0].data(), f[0].size());

            try {
                rec.flag = static_cast<std::uint16_t>(std::stoul(std::string(f[1])));
            } catch (...) {
                throw std::runtime_error("invalid SAM FLAG: " + std::string(f[1]));
            }

            rec.rname.assign(f[2].data(), f[2].size());

            try {
                rec.pos = static_cast<std::uint32_t>(std::stoul(std::string(f[3])));
            } catch (...) {
                throw std::runtime_error("invalid SAM POS: " + std::string(f[3]));
            }

            try {
                const unsigned long mapq_val = std::stoul(std::string(f[4]));
                rec.mapq = static_cast<std::uint8_t>(mapq_val > 255 ? 255 : mapq_val);
            } catch (...) {
                throw std::runtime_error("invalid SAM MAPQ: " + std::string(f[4]));
            }

            rec.cigar.assign(f[5].data(), f[5].size());
            rec.rnext.assign(f[6].data(), f[6].size());

            try {
                rec.pnext = static_cast<std::uint32_t>(std::stoul(std::string(f[7])));
            } catch (...) {
                throw std::runtime_error("invalid SAM PNEXT: " + std::string(f[7]));
            }

            try {
                rec.tlen = static_cast<std::int32_t>(std::stol(std::string(f[8])));
            } catch (...) {
                throw std::runtime_error("invalid SAM TLEN: " + std::string(f[8]));
            }

            rec.seq.assign(f[9].data(), f[9].size());
            rec.qual.assign(f[10].data(), f[10].size());

            if (!opt_fields.empty()) {
                rec.opt.assign(opt_fields.data(), opt_fields.size());
            } else {
                rec.opt.clear();
            }

            return true;
        }

        return false;
    }

    // SAM 转 FASTA
    void convertSamToFasta(const FilePath& sam_path, const FilePath& fasta_path, std::size_t line_width)
    {
        SamReader reader(sam_path);
        SeqWriter writer(fasta_path, line_width);

        SamRecord sam_rec;
        std::size_t count = 0;

        while (reader.next(sam_rec)) {
            const SeqRecord fasta_rec = samRecordToSeqRecord(sam_rec, /*keep_qual=*/false);
            writer.writeFasta(fasta_rec);
            ++count;
        }

        writer.flush();

#ifdef _DEBUG
        spdlog::debug("convertSamToFasta: converted {} records from {} to {}", count, sam_path.string(), fasta_path.string());
#endif
    }

}  // namespace seq_io

