#include "preprocess.h"
#include <chrono>
// 新增头文件，用于调用外部命令、构造字符串和检查文件存在性
#include <cstdlib>
#include <sstream>
#include <filesystem>

// 本文件包含输入 FASTA 的预处理逻辑：
// - 将输入文件（可以是本地路径或 URL）获取到工作目录的 data/raw 下；
// - 对序列做简单清洗（大写化，非 AGCTU 字符替换为 N）；
// - 将清洗后的序列写入 data/clean 下；
// - 维护一个 Top-K 选择器，选择长度最长的前 K 条序列（稳定保留较早出现的序列）；
// - 最终把清洗后的数据与选中的 consensus 输出到指定文件（在此文件中只负责选择与写出接口，具体写出由调用者处理）。

// 解释：此处使用的 FilePath 为 file_io::FilePath（即 std::filesystem::path 的别名），
// seq_io 命名空间封装了读取/写入 FASTA 的细节（openKseqReader / SeqWriter / SeqRecord 等）。

uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n) {
    // 参数说明：
    // - input_path: 输入 FASTA 的路径或 URL（字符串）。
    // - workdir: 工作目录路径（应当已经准备好或由调用方保证），在此目录下会创建 data/raw 和 data/clean 等子目录。
    // - cons_n: 要保留的最长序列数量（Top-K）。

    // 计时用于日志记录，帮助性能分析
    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Preprocessing input FASTA file: {}", input_path);
    spdlog::info("Working directory: {}", workdir);

    // ---------- 目录准备 ----------
    // 1) 确保工作目录下存在 data 文件夹。
    //    该目录用于存放原始与清洗后的数据：data/raw 和 data/clean。
    FilePath data_dir = FilePath(workdir) / WORKDIR_DATA;
    file_io::ensureDirectoryExists(data_dir);
    spdlog::info("Ensured data directory exists: {}", data_dir.string());

    // 2) 在 data 下创建 raw_data 文件夹，用于保存原始（未清洗）输入。
    FilePath raw_data_dir = data_dir / DATA_RAW;
    file_io::ensureDirectoryExists(raw_data_dir);
    spdlog::info("Ensured raw data directory exists: {}", raw_data_dir.string());

    // 3) 在 data 下创建 clean_data 文件夹，用于保存清洗后的输出。
    FilePath clean_data_dir = data_dir / DATA_CLEAN;
    file_io::ensureDirectoryExists(clean_data_dir);
    spdlog::info("Ensured clean data directory exists: {}", clean_data_dir.string());

    // ========================================================================
    // 获取输入文件（支持本地/远程，性能优化）
    // ========================================================================
    // 核心优化：
    // - 远程文件（URL）：下载到 raw_data 目录（需要本地缓存）
    // - 本地文件：直接读取原始路径（避免不必要的复制开销）
    //
    // 性能影响：
    // - 本地大文件（例如 10GB FASTA）：从数分钟复制时间减少到 0 秒
    // - 远程文件：保持原有逻辑（必须下载到本地）
    // - 内存占用：无影响（都是流式读取）
    // ========================================================================
    FilePath input_file = FilePath(input_path);
    FilePath actual_input_file;  // 实际读取的文件路径

    if (file_io::isUrl(input_file)) {
        // 远程文件：下载到 raw_data 目录
        FilePath raw_dest_file = raw_data_dir / input_file.filename();

        spdlog::info("Detected remote URL, downloading to: {} -> {}",
                     input_file.string(), raw_dest_file.string());

        // fetchFile 内部会调用 downloadFile
        file_io::fetchFile(input_file, raw_dest_file);

        spdlog::info("Download completed: {}", raw_dest_file.string());
        actual_input_file = raw_dest_file;
    } else {
        // 本地文件：直接使用原始路径，不复制
        spdlog::info("Detected local file, reading directly from: {}", input_file.string());

        // 验证文件存在性与可读性（失败时抛出异常）
        file_io::requireRegularFile(input_file, "input file");

        actual_input_file = input_file;

        spdlog::info("Local file verified, no copy needed (performance optimization)");
    }

    // ---------- 准备输出文件名 ----------
    // 5) 打开 raw 文件并逐条读取；对每条序列进行清洗（cleanSequence），写入 clean_data
    //    同时维护 TopKLongestSelector，选择最长的 cons_n 条序列（用于之后的 consensus 生成）。
    // handle input filenames like `sample.fasta.gz` -> `sample.fasta`
    FilePath in_fname = input_file.filename();
    std::string in_name = in_fname.string();
    const std::string comp_ext = ".gz";
    if (in_name.size() > comp_ext.size() &&
        in_name.compare(in_name.size() - comp_ext.size(), comp_ext.size(), comp_ext) == 0) {
        in_name.resize(in_name.size() - comp_ext.size());
        spdlog::info("Detected compressed input; using output name: {}", in_name);
    }
    FilePath clean_dest_file = clean_data_dir / FilePath(in_name);
    FilePath consensus_file = clean_data_dir / CLEAN_CONS_UNALIGNED;
    spdlog::info("Clean output: {} ; Consensus output: {}", clean_dest_file.string(), consensus_file.string());

    // ---------- 打开 reader/writer 与 TopK 选择器 ----------
    // seq_io::openKseqReader 返回一个 reader 指针（抽象），用于逐条读取序列；
    // seq_io::SeqWriter 用于把清洗后的序列写入到目标文件。
    //
    // 说明：actual_input_file 可能是：
    // 1. 远程文件：raw_data 目录下的下载文件
    // 2. 本地文件：用户提供的原始路径（无需复制，性能优化）
    auto reader = seq_io::openKseqReader(actual_input_file);
    seq_io::SeqWriter clean_writer(clean_dest_file);
    TopKLongestSelector selector(cons_n);

    // 处理循环：读取 -> 清洗 -> 写出 -> 交给 TopK 选择器
    seq_io::SeqRecord rec;
    std::size_t total_records = 0;

    // 使用 ProgressBar 类替代手动进度条实现，减少代码冗余
    ProgressBar progress("preprocess");

    // 重要说明（性能与正确性）：
    // - 该循环为预处理的热路径，若输入很大（成千上万 / 几百万条序列），要关注 IO 与内存占用。
    // - 性能优化点：
    //    * 使用 seq_io 的 KseqReader（基于 fread/gzread）并配合大缓冲可以显著提升读取吞吐；
    //    * SeqWriter::write 会把一条记录的 header 与折行后的序列缓存在临时字符串中并一次性写出，
    //      避免逐字符写入带来的系统调用开销；这对写大文件非常重要；
    //    * TopKLongestSelector 应该实现为维护一个大小为 K 的最小堆，插入/替换成本为 O(log K)，适合 K 远小于记录总数的场景；
    // - 内存权衡：TopK 的实现会保留 K 条完整记录（占用内存 O(K * avg_len)），若 K 很大需注意内存使用。
    while (reader->next(rec)) {
        ++total_records;
        // 对序列进行规范化清洗：例如把字母转为大写，非 AGCTU 替换为 N（具体实现由 seq_io::cleanSequence 提供）。
        // cleanSequence 就地修改 rec.seq，尽量避免重复复制以节省内存带宽。
        seq_io::cleanSequence(rec.seq);

        // 将清洗后的记录写入 clean_data 文件夹
        // 注：FastaWriter::write 已经在内部做了拼接与一次性写出的优化，性能友好。
        clean_writer.write(rec);

        // 将当前记录交给 TopK 选择器进行考虑（内部维护堆以保证 O(log K) 的替换成本）
        // 注意：selector.consider 应复制或接管必要的字段（例如 id/seq），以免后续 rec 被复用/覆盖导致数据错误。
        selector.consider(rec);

        // 使用 ProgressBar::tick() 更新进度，内部自动判断是否需要刷新
        progress.tick();
    }
    progress.done();  // 强制输出最终状态并换行

    // ---------- 将 TopK 结果写成共识输入文件 ----------
    // 说明：takeSortedDesc 返回按长度降序排序的记录列表（一般用于选取最长的 N 条序列作为共识计算输入）
    seq_io::SeqWriter cons_writer(consensus_file);
    auto cons_seqs = selector.takeSortedDesc();

    // 将选出的序列写到 consensus 文件；注意保持一致的换行宽度等格式规则，SeqWriter 负责这些细节。
    for (const auto& cons_rec : cons_seqs) {
        // 这里写出的 cons_rec 应该是一个深拷贝的 SeqRecord（由 selector 返回以保证安全），如果不是需要在 selector 中做拷贝。
        cons_writer.write(cons_rec);
    }

    // ---------- 统计与返回值 ----------
    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    spdlog::info("Preprocessing completed. Total records processed: {}. Selected top {} sequences: {}. Elapsed: {:.2f} s",
                 total_records, cons_n, cons_seqs.size(), elapsed_s);
    // 将 size_t 转为项目级别的 uint_t（在 config.hpp 中定义）；防止溢出则截断到 U_MAX
    if (total_records > static_cast<std::size_t>(U_MAX)) {
        spdlog::warn("Processed records ({}) exceed U_MAX ({}); truncating to U_MAX", total_records, U_MAX);
        return static_cast<uint_t>(U_MAX);
    }

    return static_cast<uint_t>(total_records);
}


void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, int threads)
{

    // 检查输入文件是否存在
    if (!std::filesystem::exists(input_file)) {
        spdlog::warn("Consensus unaligned file not found: {}", input_file.string());
        return;
    }

    // 记录开始时间
    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Starting consensus alignment");
    spdlog::info("  input : {}", input_file.string());
    spdlog::info("  output: {}", output_file.string());
    spdlog::info("  tool  : {}", msa_cmd);
    spdlog::info("  thrs  : {}", threads);

    // 尝试记录输入文件大小（若可访问）
    try {
        if (std::filesystem::exists(input_file)) {
            auto in_size = std::filesystem::file_size(input_file);
            spdlog::info("Input file size: {} bytes", in_size);
        }
    } catch (const std::exception &e) {
        spdlog::warn("Failed to stat input file {}: {}", input_file.string(), e.what());
    }

    // 组装命令。默认把 -i / -o / -t 作为参数传入，便于未来统一替换为 cmd 模块的调用。
    cmd::BuildOptions build_opt;
    const std::string cmd_str = cmd::buildCommand(msa_cmd,input_file.string(),output_file.string(),threads, build_opt);
    spdlog::info("Built MSA command (length {}): {}", cmd_str.size(), cmd_str);
    spdlog::info("MSA command (escaped): {}", cmd_str);

    // 调用外部命令（当前使用 std::system；若以后替换为项目内 cmd 接口，只需修改此处）
    try {
        const auto cmd_start = std::chrono::steady_clock::now();
        int rc = cmd::runCommand(cmd_str);
        const auto cmd_end = std::chrono::steady_clock::now();
        const double cmd_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(cmd_end - cmd_start).count();

        if (rc != 0) {
            spdlog::error("MSA command failed (exit code {}): {}", rc, cmd_str);
        } else {
            spdlog::info("MSA command exited with code 0 (success). Elapsed: {:.3f} s", cmd_elapsed);
        }

        // 检查输出文件
        try {
            if (std::filesystem::exists(output_file)) {
                auto out_size = std::filesystem::file_size(output_file);
                spdlog::info("Aligned consensus output exists: {} ({} bytes)", output_file.string(), out_size);
                if (out_size == 0) {
                    spdlog::warn("Aligned consensus output is empty: {}", input_file.string());
                }
            } else {
                spdlog::warn("Aligned consensus output not found after running MSA command: {}", input_file.string());
            }
        } catch (const std::exception &e) {
            spdlog::warn("Failed to stat output file {}: {}", input_file.string(), e.what());
        }

    } catch (const std::exception &e) {
        spdlog::error("Exception while running MSA command for {}: {}", input_file.string(), e.what());
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    spdlog::info("Finished consensus alignment. Total elapsed: {:.3f} s", elapsed_s);
}

// ==============================================================
// validateRefAlignedConsistency
//
// 说明：验证两个 FASTA 文件的序列一致性（删除 gap 后）。
// 使用 hash map 支持序列顺序不同的情况。
//
// 实现步骤：
// 1) 读取 ref_fasta，按 ID 建立 unordered_map（key=ID, value=去gap后的序列）
// 2) 逐条读取 ref_aligned，从 map 中查找匹配的序列并比较
// 3) 最后检查 ref_fasta 的序列是否全部被消耗
//
// ==============================================================
void validateRefAlignedConsistency(const FilePath& ref_fasta, const FilePath& ref_aligned)
{
    spdlog::info("Validating reference and aligned reference consistency");
    spdlog::info("  -r/--ref: {}", ref_fasta.string());
    spdlog::info("  --ref-align: {}", ref_aligned.string());

    auto toUpperBase = [](char c) -> char {
        return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    };

    // 第一步：读取 ref_fasta
    // 要求：ref_fasta 不能有 gap，并且统一转成大写保存
    std::unordered_map<std::string, std::string> ref_seqs;

    {
        seq_io::KseqReader reader(ref_fasta);
        seq_io::SeqRecord rec;

        while (reader.next(rec)) {
            if (ref_seqs.find(rec.id) != ref_seqs.end()) {
                std::ostringstream oss;
                oss << "Duplicate sequence ID in -r/--ref: '" << rec.id << "'";
                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }

            std::string seq;
            seq.reserve(rec.seq.size());

            for (std::size_t i = 0; i < rec.seq.size(); ++i) {
                char c = rec.seq[i];

                if (c == '-') {
                    std::ostringstream oss;
                    oss << "Gap character found in -r/--ref for ID '" << rec.id
                        << "' at position " << (i + 1)
                        << ". Reference FASTA must not contain gaps.";
                    spdlog::error(oss.str());
                    throw std::runtime_error(oss.str());
                }

                seq += toUpperBase(c);
            }

            ref_seqs.emplace(rec.id, std::move(seq));
        }
    }

    spdlog::info("Loaded {} sequences from -r/--ref", ref_seqs.size());

    // 第二步：读取 ref_aligned
    // 要求：
    // 1) ID 必须存在于 ref_fasta
    // 2) ref_aligned 所有序列的对齐长度必须一致
    // 3) 删除 gap 后，与 ref_fasta 对应序列一致，忽略大小写
    std::unordered_set<std::string> matched;
    std::optional<std::size_t> expected_aligned_len;

    {
        seq_io::KseqReader reader(ref_aligned);
        seq_io::SeqRecord rec;
        std::size_t seq_idx = 0;

        while (reader.next(rec)) {
            ++seq_idx;

            if (matched.find(rec.id) != matched.end()) {
                std::ostringstream oss;
                oss << "Duplicate sequence ID in --ref-align: '" << rec.id << "'";
                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }

            // 检查 aligned FASTA 中所有序列长度是否一致
            if (!expected_aligned_len.has_value()) {
                expected_aligned_len = rec.seq.size();
            } else if (rec.seq.size() != expected_aligned_len.value()) {
                std::ostringstream oss;
                oss << "Aligned sequence length mismatch in --ref-align for ID '" << rec.id
                    << "': expected aligned length " << expected_aligned_len.value()
                    << ", but got " << rec.seq.size()
                    << ". All sequences in --ref-align must have the same aligned length.";
                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }

            auto it = ref_seqs.find(rec.id);
            if (it == ref_seqs.end()) {
                std::ostringstream oss;
                oss << "Sequence ID in --ref-align not found in -r/--ref: '" << rec.id << "'";
                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }

            const std::string& ref_seq = it->second;

            // 删除 gap，并记录 ungapped position 对应的 aligned position
            std::string align_ungapped;
            std::vector<std::size_t> aligned_pos_of_ungapped_base;

            align_ungapped.reserve(rec.seq.size());
            aligned_pos_of_ungapped_base.reserve(rec.seq.size());

            for (std::size_t i = 0; i < rec.seq.size(); ++i) {
                char c = rec.seq[i];

                if (c == '-') {
                    continue;
                }

                align_ungapped += toUpperBase(c);
                aligned_pos_of_ungapped_base.push_back(i + 1);  // 1-based aligned position
            }

            // 先检查长度
            if (ref_seq.size() != align_ungapped.size()) {
                std::ostringstream oss;
                oss << "Sequence length mismatch after removing gaps for ID '" << rec.id << "': "
                    << "-r has " << ref_seq.size() << " bases, "
                    << "--ref-align has " << align_ungapped.size() << " bases.";

                const std::size_t min_len = std::min(ref_seq.size(), align_ungapped.size());
                if (min_len < ref_seq.size()) {
                    oss << " First extra base in -r is at ungapped position " << (min_len + 1)
                        << ": -r='" << ref_seq[min_len] << "', --ref-align=<end>.";
                } else if (min_len < align_ungapped.size()) {
                    oss << " First extra base in --ref-align is at ungapped position " << (min_len + 1)
                        << ", aligned position " << aligned_pos_of_ungapped_base[min_len]
                        << ": -r=<end>, --ref-align='" << align_ungapped[min_len] << "'.";
                }

                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }

            // 再检查具体碱基差异
            for (std::size_t i = 0; i < ref_seq.size(); ++i) {
                if (ref_seq[i] != align_ungapped[i]) {
                    std::ostringstream oss;
                    oss << "Sequence content mismatch after removing gaps for ID '" << rec.id << "' "
                        << "at ungapped position " << (i + 1)
                        << ", aligned position " << aligned_pos_of_ungapped_base[i]
                        << ": -r='" << ref_seq[i]
                        << "', --ref-align='" << align_ungapped[i] << "'.";

                    spdlog::error(oss.str());
                    throw std::runtime_error(oss.str());
                }
            }

            matched.insert(rec.id);
        }

        spdlog::info("Validated {} sequences from --ref-align", seq_idx);
    }

    // 第三步：检查 ref_fasta 中所有序列是否都在 ref_aligned 中出现
    if (matched.size() != ref_seqs.size()) {
        for (const auto& kv : ref_seqs) {
            if (matched.find(kv.first) == matched.end()) {
                std::ostringstream oss;
                oss << "Sequence ID in -r/--ref not found in --ref-align: '" << kv.first << "'";
                spdlog::error(oss.str());
                throw std::runtime_error(oss.str());
            }
        }

        std::ostringstream oss;
        oss << "Sequence count mismatch: -r has " << ref_seqs.size()
            << " sequences, but --ref-align has " << matched.size()
            << " matching sequences.";
        spdlog::error(oss.str());
        throw std::runtime_error(oss.str());
    }

    spdlog::info(
        "Successfully validated --ref-align consistency: all {} sequences matched",
        ref_seqs.size()
    );
}
