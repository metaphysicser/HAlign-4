#include <doctest/doctest.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "align.h"
#include "seed.h"  // 用于 minimizer 提取和锚点生成
#include "mash.h"  // 用于计算 Mash 相似度

// ------------------------------------------------------------------
// 辅助函数：生成随机 DNA 序列
// ------------------------------------------------------------------
static std::string generateRandomDNA(size_t length, unsigned seed = 42) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += bases[dist(rng)];
    }
    return seq;
}

// ------------------------------------------------------------------
// 辅助函数：在序列中引入随机突变（SNP + Indel）
// ------------------------------------------------------------------
static std::string mutateSequence(const std::string& ref,
                                  double snp_rate = 0.01,     // 1% SNP
                                  double indel_rate = 0.005,  // 0.5% indel
                                  unsigned seed = 43) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> base_dist(0, 3);
    std::uniform_int_distribution<int> indel_len(1, 5);  // indel 长度 1-5bp

    std::string mutated;
    mutated.reserve(ref.size() * 1.1);  // 预留空间

    for (size_t i = 0; i < ref.size(); ++i) {
        double p = prob(rng);

        if (p < indel_rate) {
            // 随机插入或删除
            if (prob(rng) < 0.5) {
                // 插入
                int len = indel_len(rng);
                for (int j = 0; j < len; ++j) {
                    mutated += bases[base_dist(rng)];
                }
            } else {
                // 删除：跳过当前碱基
                continue;
            }
        } else if (p < indel_rate + snp_rate) {
            // SNP：替换为不同的碱基
            char original = ref[i];
            char replacement;
            do {
                replacement = bases[base_dist(rng)];
            } while (replacement == original);
            mutated += replacement;
        } else {
            // 保持不变
            mutated += ref[i];
        }
    }

    return mutated;
}

// ------------------------------------------------------------------
// 辅助函数：生成包含结构变异的序列
// ------------------------------------------------------------------
// 支持的 SV 类型：
// - 大片段插入 (INS)
// - 大片段删除 (DEL)
// - 倒位 (INV)
// - 串联重复 (DUP)
static std::string generateSVSequence(const std::string& ref,
                                     const std::string& sv_type,
                                     size_t sv_pos,
                                     size_t sv_size,
                                     unsigned seed = 44) {
    std::string result;
    result.reserve(ref.size() + sv_size);

    if (sv_type == "INS") {
        // 插入：在指定位置插入随机序列
        result = ref.substr(0, sv_pos);
        result += generateRandomDNA(sv_size, seed);
        result += ref.substr(sv_pos);
    }
    else if (sv_type == "DEL") {
        // 删除：删除指定位置开始的片段
        result = ref.substr(0, sv_pos);
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else if (sv_type == "INV") {
        // 倒位：反转指定区域
        result = ref.substr(0, sv_pos);
        std::string inv_region = ref.substr(sv_pos, sv_size);
        std::reverse(inv_region.begin(), inv_region.end());
        // 反向互补
        for (char& c : inv_region) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
            }
        }
        result += inv_region;
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else if (sv_type == "DUP") {
        // 串联重复：复制指定区域并插入到后面
        result = ref.substr(0, sv_pos + sv_size);
        result += ref.substr(sv_pos, sv_size);  // 重复一次
        if (sv_pos + sv_size < ref.size()) {
            result += ref.substr(sv_pos + sv_size);
        }
    }
    else {
        // 未知类型，返回原序列
        result = ref;
    }

    return result;
}

// ------------------------------------------------------------------
// 辅助函数：生成包含多个 SV 的复杂序列
// ------------------------------------------------------------------
struct SVEvent {
    std::string type;  // INS, DEL, INV, DUP
    size_t pos;        // 位置
    size_t size;       // 大小
};

static std::string generateComplexSVSequence(const std::string& ref,
                                            const std::vector<SVEvent>& events,
                                            unsigned seed = 45) {
    std::string result = ref;

    // 从后往前应用 SV，避免位置偏移问题
    std::vector<SVEvent> sorted_events = events;
    std::sort(sorted_events.begin(), sorted_events.end(),
              [](const SVEvent& a, const SVEvent& b) { return a.pos > b.pos; });

    for (const auto& sv : sorted_events) {
        result = generateSVSequence(result, sv.type, sv.pos, sv.size, seed++);
    }

    return result;
}

// ------------------------------------------------------------------
// 辅助函数：生成真实的锚点（基于 minimizer 匹配）
// ------------------------------------------------------------------
// 说明：
// 之前的测试用例使用固定位置锚点（假设 ref/query 在相同位置对齐），
// 但经过 mutateSequence 后，由于 indel，这些锚点完全不准确。
//
// 本函数使用 minimizer 提取真实的共享 k-mer，生成准确的锚点：
// 1) 从 ref 和 query 中分别提取 minimizer
// 2) 使用 collect_anchors 找到匹配的 minimizer hits
// 3) 返回真实的锚点列表
//
// 参数：
// @param ref - 参考序列
// @param query - 查询序列
// @param k - k-mer 大小（默认 15）
// @param w - 窗口大小（默认 10）
// @return 真实的锚点列表
// ------------------------------------------------------------------
static anchor::Anchors generateRealAnchors(const std::string& ref,
                                           const std::string& query,
                                           std::size_t k = 15,
                                           std::size_t w = 10) {
    // 1) 提取 ref 的 minimizer
    minimizer::MinimizerHits ref_hits = minimizer::extractMinimizer(ref, k, w, false);

    // 2) 提取 query 的 minimizer
    minimizer::MinimizerHits qry_hits = minimizer::extractMinimizer(query, k, w, false);

    // 3) 收集锚点（找到匹配的 minimizer）
    anchor::Anchors anchors = minimizer::collect_anchors(ref_hits, qry_hits);

    return anchors;
}

// ------------------------------------------------------------------
// 辅助函数：CIGAR 转字符串（便于调试）
// ------------------------------------------------------------------
static std::string cigarToString(const cigar::Cigar_t& cigar) {
    std::string result;
    for (auto unit : cigar) {
        char op;
        uint32_t len;
        cigar::intToCigar(unit, op, len);
        result += std::to_string(len) + op;
    }
    return result;
}

// ------------------------------------------------------------------
// 辅助函数：写一个最小 FASTA 文件，供 CLI 解析与文件存在性校验使用。
// 说明：这里只需要“可被 ExistingFile 接受”的真实文件，不需要复杂内容。
// ------------------------------------------------------------------
static std::filesystem::path writeMiniFasta(const std::filesystem::path& p, const std::string& name) {
    std::ofstream ofs(p, std::ios::binary);
    REQUIRE_MESSAGE(ofs.good(), "cannot create fasta: " << p.string());
    ofs << ">" << name << "\n";
    ofs << "ACGTACGT\n";
    ofs.flush();
    return p;
}

// ------------------------------------------------------------------
// 性能测试辅助：计时器
// ------------------------------------------------------------------
struct Timer {
    std::chrono::high_resolution_clock::time_point start;

    Timer() : start(std::chrono::high_resolution_clock::now()) {}

    double elapsedMs() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
};

// ------------------------------------------------------------------
// 测试套件：正确性测试
// ------------------------------------------------------------------
TEST_SUITE("align") {

    TEST_CASE("ProfileMatrix - build from consensus counts keeps MSA columns") {
        consensus::ConsensusJson cj;
        cj.num_seqs = 3;
        cj.aln_len = 2;
        cj.counts.resize(2);
        cj.counts[0].a = 2;
        cj.counts[0].dash = 1;
        cj.counts[1].c = 1;
        cj.counts[1].u = 1;
        cj.counts[1].n = 1;

        align::ProfileMatrix profile = align::ProfileMatrix::fromConsensusCounts(cj);

        CHECK(profile.len == 2);
        CHECK(profile.dim == 5);
        CHECK(profile.depth == 3);
        CHECK(profile.prof[0] == 2);
        CHECK(profile.prof[1] == 0);
        CHECK(profile.prof[5 + 1] == 1);
        CHECK(profile.prof[5 + 3] == 1);
        CHECK(profile.prof[5 + 4] == 1);
    }

    TEST_CASE("globalAlignKSW2 - 精确匹配") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignKSW2(seq, seq);

        // 应该是完全匹配，CIGAR 应该是 12M
        REQUIRE(cigar.size() >= 1);
        char op;
        uint32_t len;
        cigar::intToCigar(cigar[0], op, len);
        CHECK(op == 'M');
        CHECK(len == 12);
    }

    TEST_CASE("globalAlignKSW2 - 单个错配") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTACCGTACGT";  // 第6位 T->C 错配

        auto cigar = align::globalAlignKSW2(ref, query);
        REQUIRE(cigar.size() > 0);

        // 验证 CIGAR 不为空
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignKSW2 - 单个插入") {
        std::string ref   = "ACGTACGTACGT";
        std::string query = "ACGTAACGTACGT";  // 在第5位后插入 A

        auto cigar = align::globalAlignKSW2(ref, query);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR: ", cigar_str);

        // 应该包含插入操作（I）
        bool has_insertion = cigar_str.find('I') != std::string::npos;
        CHECK(has_insertion);
    }

    TEST_CASE("extendAlignKSW2 - 基本延伸") {
        std::string ref = generateRandomDNA(1000, 100);
        std::string query = ref.substr(100, 500);  // 提取中间片段

        auto cigar = align::extendAlignKSW2(ref, query, 200);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("Extend CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignWFA2 - 精确匹配") {
        std::string seq = "ACGTACGTACGT";
        auto cigar = align::globalAlignWFA2(seq, seq);

        REQUIRE(cigar.size() >= 1);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("WFA2 CIGAR: ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("空序列边界测试") {
        std::string empty = "";
        std::string seq = "ACGT";

        // 测试空序列不崩溃（具体行为依赖于实现）
        CHECK_NOTHROW(align::globalAlignKSW2(empty, seq));
        CHECK_NOTHROW(align::globalAlignKSW2(seq, empty));
    }

    TEST_CASE("高相似度序列 - 99% 相似") {
        // 生成一个 1000bp 的序列，引入 1% 的错配
        std::string ref = generateRandomDNA(1000, 12345);
        std::string query = mutateSequence(ref, 0.01, 0.0, 12346);  // 只有 SNP，无 indel

        // 三种方法都应该能正确比对
        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_extend = align::extendAlignKSW2(ref, query, 200);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // 验证 CIGAR 都不为空
        CHECK(cigar_ksw2.size() > 0);
        CHECK(cigar_extend.size() > 0);
        CHECK(cigar_wfa2.size() > 0);

        MESSAGE("KSW2 CIGAR: ", cigarToString(cigar_ksw2));
        MESSAGE("Extend CIGAR: ", cigarToString(cigar_extend));
        MESSAGE("WFA2 CIGAR: ", cigarToString(cigar_wfa2));
    }

    TEST_CASE("高相似度序列 - 98% 相似（含 indel）") {
        // 生成序列，引入 1% SNP + 1% indel
        std::string ref = generateRandomDNA(500, 54321);
        std::string query = mutateSequence(ref, 0.01, 0.01, 54322);

        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);
        auto cigar_wfa2 = align::globalAlignWFA2(ref, query);

        // 验证 CIGAR 包含不同类型的操作
        std::string cigar_str_ksw2 = cigarToString(cigar_ksw2);
        std::string cigar_str_wfa2 = cigarToString(cigar_wfa2);

        MESSAGE("KSW2 CIGAR (with indels): ", cigar_str_ksw2);
        MESSAGE("WFA2 CIGAR (with indels): ", cigar_str_wfa2);

        CHECK(!cigar_str_ksw2.empty());
        CHECK(!cigar_str_wfa2.empty());

        // 应该包含匹配操作
        bool has_match_ksw2 = cigar_str_ksw2.find('M') != std::string::npos;
        bool has_match_wfa2 = cigar_str_wfa2.find('M') != std::string::npos;
        CHECK(has_match_ksw2);
        CHECK(has_match_wfa2);
    }

    TEST_CASE("极高相似度序列 - 99.9% 相似") {
        // 模拟测序错误：仅 0.1% 的错误率
        std::string ref = generateRandomDNA(10000, 99999);
        std::string query = mutateSequence(ref, 0.0005, 0.0005, 100000);

        // 在如此高的相似度下，所有方法都应该快速完成
        Timer timer;
        auto cigar = align::globalAlignKSW2(ref, query);
        double elapsed = timer.elapsedMs();

        CHECK(cigar.size() > 0);
        MESSAGE("极高相似度 10k 序列比对耗时: ", elapsed, " ms");

        // 对于 99.9% 相似的 10k 序列，应该在合理时间内完成（<100ms）
        CHECK(elapsed < 100.0);
    }

    // ------------------------------------------------------------------
    // 测试：cigarToString 和 stringToCigar 互逆
    // ------------------------------------------------------------------
    TEST_CASE("cigar::cigarToString and stringToCigar - 互逆操作") {
        // 测试1：标准 CIGAR 字符串
        SUBCASE("标准 CIGAR") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 100));
            original.push_back(cigar::cigarToInt('I', 5));
            original.push_back(cigar::cigarToInt('M', 95));
            original.push_back(cigar::cigarToInt('D', 3));
            original.push_back(cigar::cigarToInt('M', 50));

            // 转换为字符串
            std::string cigar_str = cigar::cigarToString(original);
            CHECK(cigar_str == "100M5I95M3D50M");

            // 再转回 Cigar_t
            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);

            // 验证互逆
            REQUIRE(roundtrip.size() == original.size());
            for (size_t i = 0; i < original.size(); ++i) {
                CHECK(roundtrip[i] == original[i]);
            }
        }

        // 测试2：所有 CIGAR 操作符
        SUBCASE("所有操作符") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 10));
            original.push_back(cigar::cigarToInt('I', 2));
            original.push_back(cigar::cigarToInt('D', 3));
            original.push_back(cigar::cigarToInt('N', 100));
            original.push_back(cigar::cigarToInt('S', 5));
            original.push_back(cigar::cigarToInt('H', 10));
            original.push_back(cigar::cigarToInt('P', 1));
            original.push_back(cigar::cigarToInt('=', 20));
            original.push_back(cigar::cigarToInt('X', 3));

            std::string cigar_str = cigar::cigarToString(original);
            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);

            REQUIRE(roundtrip.size() == original.size());
            for (size_t i = 0; i < original.size(); ++i) {
                CHECK(roundtrip[i] == original[i]);
            }
        }

        // 测试3：特殊值 "*"
        SUBCASE("特殊值 *") {
            cigar::Cigar_t empty_cigar = cigar::stringToCigar("*");
            CHECK(empty_cigar.empty());
        }

        // 测试4：空字符串
        SUBCASE("空字符串") {
            cigar::Cigar_t empty_cigar = cigar::stringToCigar("");
            CHECK(empty_cigar.empty());
        }

        // 测试5：大数字长度
        SUBCASE("大数字长度") {
            cigar::Cigar_t original;
            original.push_back(cigar::cigarToInt('M', 999999));
            original.push_back(cigar::cigarToInt('D', 123456));

            std::string cigar_str = cigar::cigarToString(original);
            CHECK(cigar_str == "999999M123456D");

            cigar::Cigar_t roundtrip = cigar::stringToCigar(cigar_str);
            REQUIRE(roundtrip.size() == original.size());
            CHECK(roundtrip[0] == original[0]);
            CHECK(roundtrip[1] == original[1]);
        }
    }

    // ------------------------------------------------------------------
    // 测试：stringToCigar 错误处理
    // ------------------------------------------------------------------
    TEST_CASE("cigar::stringToCigar - 错误处理") {
        // 测试1：操作符前没有数字
        SUBCASE("操作符前没有数字") {
            CHECK_THROWS_AS(cigar::stringToCigar("M10"), std::runtime_error);
        }

        // 测试2：未知操作符
        SUBCASE("未知操作符") {
            CHECK_THROWS_AS(cigar::stringToCigar("10Q"), std::runtime_error);
        }

        // 测试3：字符串结尾有数字但没有操作符
        SUBCASE("结尾有数字无操作符") {
            CHECK_THROWS_AS(cigar::stringToCigar("10M5"), std::runtime_error);
        }

        // 测试4：长度为 0
        SUBCASE("长度为 0") {
            CHECK_THROWS_AS(cigar::stringToCigar("0M"), std::runtime_error);
        }
    }

    // ------------------------------------------------------------------
    // 测试：stringToCigar 容错性
    // ------------------------------------------------------------------
    TEST_CASE("cigar::stringToCigar - 容错性") {
        // 测试1：包含空白字符
        SUBCASE("包含空白字符") {
            std::string cigar_with_spaces = " 10M 5I  3D ";
            cigar::Cigar_t result = cigar::stringToCigar(cigar_with_spaces);

            REQUIRE(result.size() == 3);
            char op;
            uint32_t len;

            cigar::intToCigar(result[0], op, len);
            CHECK(op == 'M');
            CHECK(len == 10);

            cigar::intToCigar(result[1], op, len);
            CHECK(op == 'I');
            CHECK(len == 5);

            cigar::intToCigar(result[2], op, len);
            CHECK(op == 'D');
            CHECK(len == 3);
        }
    }

    TEST_CASE("CLI setup - parse --reference-msa") {
        namespace fs = std::filesystem;

        const fs::path dir = fs::current_path() / "halign4_tests_cli_ref_align";
        std::error_code ec;
        fs::remove_all(dir, ec);
        fs::create_directories(dir, ec);
        REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");

        const fs::path input = writeMiniFasta(dir / "input.fasta", "input");
        const fs::path ref = writeMiniFasta(dir / "ref.fasta", "ref");
        const fs::path ref_align = writeMiniFasta(dir / "ref_aligned.fasta", "ref_aligned");
        const fs::path output = dir / "out.fasta";

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);

        std::vector<std::string> args = {
            "halign4",
            "-i", input.string(),
            "-o", output.string(),
            "-r", ref.string(),
            "--reference-msa", ref_align.string()
        };
        std::vector<char*> argv;
        argv.reserve(args.size());
        for (auto& s : args) {
            argv.push_back(s.data());
        }

        REQUIRE_NOTHROW(app.parse(static_cast<int>(argv.size()), argv.data()));
        CHECK(opt.reference_path == ref.string());
        CHECK(opt.reference_msa_path == ref_align.string());

        fs::remove_all(dir, ec);
    }

    TEST_CASE("CLI setup - parse aligned reference flag") {
        namespace fs = std::filesystem;

        const fs::path dir = fs::current_path() / "halign4_tests_cli_ref_aligned_flag";
        std::error_code ec;
        fs::remove_all(dir, ec);
        fs::create_directories(dir, ec);
        REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");

        const fs::path input = writeMiniFasta(dir / "input.fasta", "input");
        const fs::path ref_align = writeMiniFasta(dir / "ref_aligned.fasta", "ref_aligned");
        const fs::path output = dir / "out.fasta";

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);

        std::vector<std::string> args = {
            "halign4",
            "-i", input.string(),
            "-o", output.string(),
            "-r", ref_align.string(),
            "-a"
        };
        std::vector<char*> argv;
        argv.reserve(args.size());
        for (auto& s : args) {
            argv.push_back(s.data());
        }

        REQUIRE_NOTHROW(app.parse(static_cast<int>(argv.size()), argv.data()));
        CHECK(opt.reference_path == ref_align.string());
        CHECK(opt.reference_is_aligned);
        CHECK(opt.reference_msa_path.empty());

        fs::remove_all(dir, ec);
    }

    TEST_CASE("CLI setup - parse profile reference options") {
        namespace fs = std::filesystem;

        const fs::path dir = fs::current_path() / "halign4_tests_cli_profile_ref_options";
        std::error_code ec;
        fs::remove_all(dir, ec);
        fs::create_directories(dir, ec);
        REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");

        const fs::path input = writeMiniFasta(dir / "input.fasta", "input");
        const fs::path output = dir / "out.fasta";

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);

        std::vector<std::string> args = {
            "halign4",
            "-i", input.string(),
            "-o", output.string(),
            "--sketch-kmer-size", "10",
            "--min-profile-references", "15",
            "--max-profile-references", "40",
            "--min-profile-reference-similarity", "0.7"
        };
        std::vector<char*> argv;
        argv.reserve(args.size());
        for (auto& s : args) {
            argv.push_back(s.data());
        }

        REQUIRE_NOTHROW(app.parse(static_cast<int>(argv.size()), argv.data()));
        CHECK(opt.sketch_kmer_size == 10);
        CHECK(opt.min_profile_references == 15);
        CHECK(opt.max_profile_references == 40);
        CHECK(opt.min_profile_reference_similarity == doctest::Approx(0.7));

        fs::remove_all(dir, ec);
    }

    TEST_CASE("CLI setup - parse merge output options") {
        namespace fs = std::filesystem;

        const fs::path dir = fs::current_path() / "halign4_tests_cli_merge_output_options";
        std::error_code ec;
        fs::remove_all(dir, ec);
        fs::create_directories(dir, ec);
        REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");

        const fs::path input = writeMiniFasta(dir / "input.fasta", "input");
        const fs::path output = dir / "out.fasta";
        const fs::path insertion_tsv = dir / "insertions.tsv";

        Options opt;
        CLI::App app{"halign4"};
        setupCli(app, opt);

        std::vector<std::string> args = {
            "halign4",
            "-i", input.string(),
            "-o", output.string(),
            "--no-reference-output",
            "--insertions-output", insertion_tsv.string(),
            "--insertion-merge", "msa"
        };
        std::vector<char*> argv;
        argv.reserve(args.size());
        for (auto& s : args) {
            argv.push_back(s.data());
        }

        REQUIRE_NOTHROW(app.parse(static_cast<int>(argv.size()), argv.data()));
        CHECK(opt.no_reference_output);
        CHECK(opt.insertions_output == insertion_tsv.string());
        CHECK(opt.insertion_merge == "msa");

        fs::remove_all(dir, ec);
    }
}

// ------------------------------------------------------------------
// 性能测试套件
// ------------------------------------------------------------------
TEST_SUITE("align_perf") {

    // ------------------------------------------------------------------
    // 性能测试：短序列（~100bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Short sequences (~100bp)") {
        constexpr int NUM_RUNS = 1000;
        constexpr size_t SEQ_LEN = 100;

        std::cout << "\n========== 短序列性能测试 (100bp, " << NUM_RUNS << " 次) ==========\n";

        // 生成测试数据
        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // 为 MM2 生成模拟锚点（每 30bp 一个）
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 15 < SEQ_LEN; pos += 30) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 15;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        // 测试 globalAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;  // 防止优化掉
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        // 测试 extendAlignKSW2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        // 测试 globalAlignWFA2
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        // 测试 globalAlignMM2（有锚点）
        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：中等长度序列（~1000bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Medium sequences (~1000bp)") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 中等序列性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // 为 MM2 生成模拟锚点（每 150bp 一个）
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 20 < SEQ_LEN; pos += 150) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 20;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：长序列（~10000bp）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Long sequences (~10000bp)") {
        constexpr int NUM_RUNS = 10;
        constexpr size_t SEQ_LEN = 10000;

        std::cout << "\n========== 长序列性能测试 (10000bp, " << NUM_RUNS << " 次) ==========\n";

        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
            std::string query = mutateSequence(ref, 0.02, 0.01, i * 2 + 1);
            test_pairs.emplace_back(ref, query);

            // 为 MM2 生成模拟锚点（每 500bp 一个）
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += 500) {
                anchor::Anchor a;
                a.hash = pos * 1000 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);
                a.span = 50;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::extendAlignKSW2(ref, query, 200);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  extendAlignKSW2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignWFA2(ref, query);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignWFA2:  " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();
            std::cout << "  globalAlignMM2:   " << std::fixed << std::setprecision(2)
                      << elapsed << " ms (" << (elapsed / NUM_RUNS) << " ms/次)\n";
        }

        std::cout << "========================================================\n\n";
    }


    // ------------------------------------------------------------------
    // 性能测试：高相似度序列（95%-99.9%，真实测序场景）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - High similarity sequences (real-world)") {
        constexpr int NUM_RUNS = 200;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 高相似度序列性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";
        std::cout << "说明：模拟真实基因组测序场景，序列相似度 >95%\n\n";

        // 测试不同相似度等级（通过控制 SNP + indel 率）
        struct SimilarityLevel {
            double snp_rate;
            double indel_rate;
            const char* desc;
            double similarity;  // 预期相似度
        };

        std::vector<SimilarityLevel> levels = {
            {0.0001, 0.0001, "极高相似度 (>99.9%)", 99.98},
            {0.001,  0.0005, "很高相似度 (~99%)",   99.0},
            {0.005,  0.002,  "高相似度 (~98%)",     98.0},
            {0.01,   0.005,  "中等相似度 (~97%)",   97.0},
            {0.03,   0.01,   "较低相似度 (~95%)",   95.0},
            {0.05,   0.02,   "低相似度 (~90%)",     90.0},
            {0.08,   0.03,   "更低相似度 (~85%)",   85.0},
            {0.12,   0.05,   "很低相似度 (~80%)",   80.0},
            {0.18,   0.07,   "极低相似度 (~70%)",   70.0}
        };

        for (const auto& level : levels) {
            std::cout << "---------- " << level.desc << " ----------\n";

            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, level.snp_rate, level.indel_rate, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // 根据相似度调整锚点生成策略
                // 高相似度：锚点密集（每 150bp）
                // 低相似度：锚点稀疏（每 300bp），因为可靠锚点更少
                size_t anchor_interval = (level.similarity >= 90.0) ? 150 : 300;
                size_t anchor_span = (level.similarity >= 90.0) ? 20 : 15;

                anchor::Anchors anchors;
                for (size_t pos = 0; pos + anchor_span < SEQ_LEN; pos += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = static_cast<uint32_t>(anchor_span);
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            // ========== Mash 相似度计算 ==========
            // 使用 Mash sketch 计算实际的 Jaccard 相似度和 ANI
            {
                constexpr std::size_t MASH_K = 21;          // k-mer 大小
                constexpr std::size_t MASH_SKETCH_SIZE = 2000; // sketch 大小

                double total_jaccard = 0.0;
                double total_ani = 0.0;
                int valid_count = 0;

                // 对前 20 对序列计算 Mash 相似度（避免计算开销过大）
                int mash_sample_size = std::min(20, NUM_RUNS);
                for (int i = 0; i < mash_sample_size; ++i) {
                    const auto& [ref, query] = test_pairs[i];

                    // 生成 sketch
                    auto sketch_ref = mash::sketchFromSequence(ref, MASH_K, MASH_SKETCH_SIZE);
                    auto sketch_query = mash::sketchFromSequence(query, MASH_K, MASH_SKETCH_SIZE);

                    // 计算 Jaccard 相似度
                    double j = mash::jaccard(sketch_ref, sketch_query);

                    // 计算 ANI（Average Nucleotide Identity）
                    double ani = mash::aniFromJaccard(j, MASH_K);

                    total_jaccard += j;
                    total_ani += ani;
                    valid_count++;
                }

                if (valid_count > 0) {
                    double avg_jaccard = total_jaccard / valid_count;
                    double avg_ani = total_ani / valid_count;

                    std::cout << "  Mash 相似度 (k=" << MASH_K << ", s=" << MASH_SKETCH_SIZE << "):\n";
                    std::cout << "    Jaccard: " << std::fixed << std::setprecision(4) << avg_jaccard
                              << "  ANI: " << std::setprecision(2) << (avg_ani * 100.0) << "%\n";
                }
            }

            // ========== 比对性能测试 ==========

            // KSW2 测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 全局比对:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            // KSW2 延伸模式测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  KSW2 延伸模式:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            // WFA2 测试
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    // 使用预期相似度（已在 SimilarityLevel 中定义）
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  WFA2 全局比对:  " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            // MM2 测试
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                double elapsed = timer.elapsedMs();
                double avg_ms = elapsed / NUM_RUNS;
                std::cout << "  MM2 锚点比对:   " << std::fixed << std::setprecision(3)
                          << avg_ms << " ms/次  (吞吐: " << (1000.0 / avg_ms)
                          << " 次/秒)\n";
            }

            std::cout << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：长度敏感性（高相似度，不同长度）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Length scaling (high similarity)") {
        constexpr int NUM_RUNS = 50;
        constexpr double SNP_RATE = 0.01;    // 1% SNP（~98% 相似度）
        constexpr double INDEL_RATE = 0.005; // 0.5% indel

        std::cout << "\n========== 长度扩展性测试 (相似度 ~98%, " << NUM_RUNS << " 次) ==========\n";

        std::vector<size_t> lengths = {100, 500, 1000, 5000, 10000};

        std::cout << std::setw(10) << "长度(bp)"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(70, '-') << "\n";

        for (size_t len : lengths) {
            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(len, i * 2);
                std::string query = mutateSequence(ref, SNP_RATE, INDEL_RATE, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // 生成锚点（根据长度调整密度）
                size_t anchor_interval = std::max(size_t(50), len / 10);  // 大约 10 个锚点
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 20 < len; pos += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 20;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0;
            double extend_time = 0.0;
            double wfa2_time = 0.0;
            double mm2_time = 0.0;

            // KSW2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // Extend
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(10) << len
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：长度差异扩展性（固定长 ref，变化短 query）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Length difference scalability (30k ref vs varying query)") {
        constexpr int NUM_RUNS = 30;
        constexpr size_t REF_LEN = 30000;  // 固定 ref 长度为 3w bp
        constexpr double SIMILARITY = 0.95; // 95% 相似度
        constexpr double SNP_RATE = 0.03;   // 3% SNP
        constexpr double INDEL_RATE = 0.02; // 2% indel

        std::cout << "\n========== 长度差异扩展性测试 (ref=" << REF_LEN << "bp, 相似度 ~95%, " << NUM_RUNS << " 次) ==========\n";
        std::cout << "说明：模拟长 ref 序列与不同长度 query 的比对性能（测试长度不对称场景）\n\n";

        std::vector<size_t> query_lengths = {100, 500, 1000, 5000, 10000, 15000, 20000, 25000};

        std::cout << std::setw(12) << "Query长度"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(72, '-') << "\n";

        for (size_t query_len : query_lengths) {
            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                // 固定生成 30k ref
                std::string ref = generateRandomDNA(REF_LEN, i * 7);

                // 从 ref 中选择一段连续区域作为 query 的模板（模拟真实比对场景）
                // 选择中间偏前的位置，确保有足够空间
                size_t start_pos = (REF_LEN - query_len) / 3;
                std::string ref_segment = ref.substr(start_pos, query_len);

                // 在该段上施加突变，得到 query（95% 相似度）
                std::string query = mutateSequence(ref_segment, SNP_RATE, INDEL_RATE, i * 7 + 1);

                test_pairs.emplace_back(ref, query);

                // 生成锚点：在 ref 和 query 的对应区域
                // 锚点间隔根据 query 长度调整
                size_t anchor_interval = std::max(size_t(100), query_len / 8);
                anchor::Anchors anchors;

                for (size_t offset = 0; offset + 20 < query_len; offset += anchor_interval) {
                    anchor::Anchor a;
                    a.hash = (start_pos + offset) * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(start_pos + offset);  // ref 上的实际位置
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(offset);              // query 上的对应位置
                    a.span = 20;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0;
            double extend_time = 0.0;
            double wfa2_time = 0.0;
            double mm2_time = 0.0;

            // KSW2 全局比对
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // KSW2 延伸模式（从锚点延伸）
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2 全局比对
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2 锚点辅助比对
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(12) << query_len
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "\n备注：\n";
        std::cout << "  - 该测试模拟真实场景：长参考序列与变长查询序列比对\n";
        std::cout << "  - Query 是从 Ref 的某个区域提取并突变得到（95% 相似度）\n";
        std::cout << "  - 锚点辅助算法（MM2/Extend）在此场景下应显著减少搜索空间\n";
        std::cout << "  - KSW2 全局比对复杂度为 O(m*n)，在 ref 很长时会有明显开销\n";
        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：低相似度序列（70%-90%，挑战性场景）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - Low similarity sequences (70%-90%)") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 1000;

        std::cout << "\n========== 低相似度序列性能测试 (1000bp, " << NUM_RUNS << " 次) ==========\n";
        std::cout << "说明：测试算法在高突变率序列上的性能\n\n";

        struct LowSimilarityLevel {
            double snp_rate;
            double indel_rate;
            const char* desc;
        };

        std::vector<LowSimilarityLevel> levels = {
            {0.05,  0.02,  "90% 相似度"},
            {0.10,  0.05,  "85% 相似度"},
            {0.15,  0.08,  "75% 相似度"},
            {0.20,  0.10,  "70% 相似度"}
        };

        std::cout << std::setw(20) << "相似度等级"
                  << std::setw(15) << "KSW2(ms)"
                  << std::setw(15) << "Extend(ms)"
                  << std::setw(15) << "WFA2(ms)"
                  << std::setw(15) << "MM2(ms)" << "\n";
        std::cout << std::string(80, '-') << "\n";

        for (const auto& level : levels) {
            // 生成测试数据
            std::vector<std::pair<std::string, std::string>> test_pairs;
            std::vector<anchor::Anchors> anchors_list;

            for (int i = 0; i < NUM_RUNS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 2);
                std::string query = mutateSequence(ref, level.snp_rate, level.indel_rate, i * 2 + 1);
                test_pairs.emplace_back(ref, query);

                // 低相似度场景：锚点更稀疏（每 300bp 一个）
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 15 < SEQ_LEN; pos += 300) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + i;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 15;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                anchors_list.push_back(anchors);
            }

            double ksw2_time = 0.0, extend_time = 0.0, wfa2_time = 0.0, mm2_time = 0.0;

            // KSW2
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    (void)cigar;
                }
                ksw2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // Extend
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    auto cigar = align::extendAlignKSW2(ref, query, 200);
                    (void)cigar;
                }
                extend_time = timer.elapsedMs() / NUM_RUNS;
            }

            // WFA2（使用估计的相似度）
            {
                Timer timer;
                for (const auto& [ref, query] : test_pairs) {
                    double similarity = 1.0 - level.snp_rate - level.indel_rate;
                    auto cigar = align::globalAlignWFA2(ref, query);
                    (void)cigar;
                }
                wfa2_time = timer.elapsedMs() / NUM_RUNS;
            }

            // MM2
            {
                Timer timer;
                for (size_t i = 0; i < test_pairs.size(); ++i) {
                    const auto& [ref, query] = test_pairs[i];
                    auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                    (void)cigar;
                }
                mm2_time = timer.elapsedMs() / NUM_RUNS;
            }

            std::cout << std::setw(20) << level.desc
                      << std::setw(15) << std::fixed << std::setprecision(3) << ksw2_time
                      << std::setw(15) << extend_time
                      << std::setw(15) << wfa2_time
                      << std::setw(15) << mm2_time << "\n";
        }

        std::cout << "\n备注：低相似度时，MM2 的锚点可能不够可靠，性能优势会减弱\n";
        std::cout << "========================================================\n\n";
    }


    // ------------------------------------------------------------------
    // 性能测试：globalAlignMM2 vs globalAlignKSW2（带锚点 vs 无锚点）
    // ------------------------------------------------------------------
    TEST_CASE("Performance - globalAlignMM2 with anchors") {
        constexpr int NUM_RUNS = 100;
        constexpr size_t SEQ_LEN = 2000;

        std::cout << "\n========== globalAlignMM2 性能测试 (2000bp, " << NUM_RUNS << " 次) ==========\n";
        std::cout << "说明：对比有锚点辅助的 MM2 比对 vs 纯全局比对\n\n";

        // 生成测试数据：高相似度序列（98%）
        std::vector<std::pair<std::string, std::string>> test_pairs;
        std::vector<anchor::Anchors> anchors_list;

        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
            std::string query = mutateSequence(ref, 0.01, 0.01, i * 3 + 1);
            test_pairs.emplace_back(ref, query);

            // 为每对序列生成模拟锚点
            // 在高相似度序列中，每隔 200bp 创建一个锚点
            anchor::Anchors anchors;
            for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += 200) {
                anchor::Anchor a;
                a.hash = pos * 100 + i;
                a.rid_ref = 0;
                a.pos_ref = static_cast<uint32_t>(pos);
                a.rid_qry = 0;
                a.pos_qry = static_cast<uint32_t>(pos);  // 假设对齐良好
                a.span = 50;
                a.is_rev = false;
                anchors.push_back(a);
            }
            anchors_list.push_back(anchors);
        }

        // 测试 globalAlignKSW2（无锚点，纯全局比对）
        double ksw2_time = 0.0;
        {
            Timer timer;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignKSW2(ref, query);
                (void)cigar;
            }
            ksw2_time = timer.elapsedMs();
        }

        // 测试 globalAlignMM2（有锚点）
        double mm2_time = 0.0;
        {
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            mm2_time = timer.elapsedMs();
        }

        // 测试 globalAlignMM2（空锚点，应退化为 KSW2）
        double mm2_empty_time = 0.0;
        {
            Timer timer;
            anchor::Anchors empty;
            for (const auto& [ref, query] : test_pairs) {
                auto cigar = align::globalAlignMM2(ref, query, empty);
                (void)cigar;
            }
            mm2_empty_time = timer.elapsedMs();
        }

        std::cout << "  globalAlignKSW2:              " << std::fixed << std::setprecision(2)
                  << ksw2_time << " ms (" << (ksw2_time / NUM_RUNS) << " ms/次)\n";
        std::cout << "  globalAlignMM2 (有锚点):      " << std::fixed << std::setprecision(2)
                  << mm2_time << " ms (" << (mm2_time / NUM_RUNS) << " ms/次)\n";
        std::cout << "  globalAlignMM2 (空锚点):      " << std::fixed << std::setprecision(2)
                  << mm2_empty_time << " ms (" << (mm2_empty_time / NUM_RUNS) << " ms/次)\n";

        double speedup = ksw2_time / mm2_time;
        std::cout << "\n  加速比 (有锚点 vs 纯KSW2):    " << std::fixed << std::setprecision(2)
                  << speedup << "x\n";

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 性能测试：不同锚点密度下的 globalAlignMM2 性能
    // ------------------------------------------------------------------
    TEST_CASE("Performance - globalAlignMM2 anchor density") {
        constexpr int NUM_RUNS = 50;
        constexpr size_t SEQ_LEN = 3000;

        std::cout << "\n========== globalAlignMM2 锚点密度测试 (3000bp, " << NUM_RUNS << " 次) ==========\n";

        // 生成测试数据
        std::vector<std::pair<std::string, std::string>> test_pairs;
        for (int i = 0; i < NUM_RUNS; ++i) {
            std::string ref = generateRandomDNA(SEQ_LEN, i * 4);
            std::string query = mutateSequence(ref, 0.01, 0.005, i * 4 + 1);
            test_pairs.emplace_back(ref, query);
        }

        // 测试不同的锚点间隔（密度）
        std::vector<size_t> anchor_intervals = {500, 300, 200, 100, 50};

        std::cout << std::setw(15) << "锚点间隔(bp)"
                  << std::setw(15) << "锚点数量"
                  << std::setw(15) << "耗时(ms)" << "\n";
        std::cout << std::string(45, '-') << "\n";

        for (size_t interval : anchor_intervals) {
            // 生成对应密度的锚点
            std::vector<anchor::Anchors> anchors_list;
            size_t avg_anchor_count = 0;

            for (size_t run = 0; run < test_pairs.size(); ++run) {
                anchor::Anchors anchors;
                for (size_t pos = 0; pos + 50 < SEQ_LEN; pos += interval) {
                    anchor::Anchor a;
                    a.hash = pos * 1000 + run;
                    a.rid_ref = 0;
                    a.pos_ref = static_cast<uint32_t>(pos);
                    a.rid_qry = 0;
                    a.pos_qry = static_cast<uint32_t>(pos);
                    a.span = 50;
                    a.is_rev = false;
                    anchors.push_back(a);
                }
                avg_anchor_count += anchors.size();
                anchors_list.push_back(anchors);
            }
            avg_anchor_count /= test_pairs.size();

            // 测试性能
            Timer timer;
            for (size_t i = 0; i < test_pairs.size(); ++i) {
                const auto& [ref, query] = test_pairs[i];
                auto cigar = align::globalAlignMM2(ref, query, anchors_list[i]);
                (void)cigar;
            }
            double elapsed = timer.elapsedMs();

            std::cout << std::setw(15) << interval
                      << std::setw(15) << avg_anchor_count
                      << std::setw(15) << std::fixed << std::setprecision(2)
                      << (elapsed / NUM_RUNS) << "\n";
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 测试：globalAlignMM2 - 基于锚点的全局比对
    // ------------------------------------------------------------------
    TEST_CASE("globalAlignMM2 - 空锚点退化为全局比对") {
        std::string ref = "ACGTACGTACGT";
        std::string query = "ACGTACCGTACGT";  // 第6位有错配

        anchor::Anchors empty_anchors;
        auto cigar = align::globalAlignMM2(ref, query, empty_anchors);

        // 应该退化为 globalAlignKSW2
        REQUIRE(cigar.size() > 0);
        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (empty anchors): ", cigar_str);
        CHECK(!cigar_str.empty());
    }

    TEST_CASE("globalAlignMM2 - 单个锚点") {
        std::string ref = "ACGTACGTACGT";
        std::string query = "ACGTACGTACGT";

        // 创建一个锚点：在位置 4，长度 4
        anchor::Anchors anchors;
        anchor::Anchor a;
        a.hash = 12345;
        a.rid_ref = 0;
        a.pos_ref = 4;
        a.rid_qry = 0;
        a.pos_qry = 4;
        a.span = 4;
        a.is_rev = false;
        anchors.push_back(a);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (single anchor): ", cigar_str);

        // 验证 CIGAR 消耗的序列长度
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - 多个锚点形成链") {
        std::string ref = "ACGTACGTACGTACGTACGT";  // 20bp
        std::string query = "ACGTACGTACGTACGTACGT"; // 完全匹配

        // 创建多个锚点
        anchor::Anchors anchors;

        // 锚点 1: pos=0, span=4
        anchor::Anchor a1;
        a1.hash = 1001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 4;
        a1.is_rev = false;
        anchors.push_back(a1);

        // 锚点 2: pos=8, span=4
        anchor::Anchor a2;
        a2.hash = 1002;
        a2.rid_ref = 0;
        a2.pos_ref = 8;
        a2.rid_qry = 0;
        a2.pos_qry = 8;
        a2.span = 4;
        a2.is_rev = false;
        anchors.push_back(a2);

        // 锚点 3: pos=16, span=4
        anchor::Anchor a3;
        a3.hash = 1003;
        a3.rid_ref = 0;
        a3.pos_ref = 16;
        a3.rid_qry = 0;
        a3.pos_qry = 16;
        a3.span = 4;
        a3.is_rev = false;
        anchors.push_back(a3);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (multiple anchors): ", cigar_str);

        // 验证完整覆盖
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - 锚点间有间隙（小间隙）") {
        std::string ref =   "AAAA----CCCC----GGGG";  // 20bp (去掉 '-' 后 12bp)
        std::string query = "AAAATTTTCCCCTTTTGGGG";  // 20bp

        // 实际序列（无 gap）
        std::string ref_actual = "AAAACCCCGGGG";  // 12bp
        std::string query_actual = "AAAATTTTCCCCTTTTGGGG";  // 20bp

        // 创建锚点：只在匹配的区域
        anchor::Anchors anchors;

        // 锚点 1: AAAA (ref: 0-3, query: 0-3)
        anchor::Anchor a1;
        a1.hash = 2001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 4;
        a1.is_rev = false;
        anchors.push_back(a1);

        // 锚点 2: CCCC (ref: 4-7, query: 8-11)
        anchor::Anchor a2;
        a2.hash = 2002;
        a2.rid_ref = 0;
        a2.pos_ref = 4;
        a2.rid_qry = 0;
        a2.pos_qry = 8;
        a2.span = 4;
        a2.is_rev = false;
        anchors.push_back(a2);

        // 锚点 3: GGGG (ref: 8-11, query: 16-19)
        anchor::Anchor a3;
        a3.hash = 2003;
        a3.rid_ref = 0;
        a3.pos_ref = 8;
        a3.rid_qry = 0;
        a3.pos_qry = 16;
        a3.span = 4;
        a3.is_rev = false;
        anchors.push_back(a3);

        auto cigar = align::globalAlignMM2(ref_actual, query_actual, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (with small gaps): ", cigar_str);

        // 验证完整覆盖
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref_actual.size());
        CHECK(qry_len == query_actual.size());

        // 应该包含插入操作（query 比 ref 长）
        bool has_insertion = cigar_str.find('I') != std::string::npos;
        CHECK(has_insertion);
    }

    TEST_CASE("globalAlignMM2 - 锚点间有大间隙（测试自适应策略）") {
        // 创建一个 1000bp 的参考序列
        std::string ref = generateRandomDNA(1000, 5000);
        std::string query = ref;  // 先完全匹配

        // 在 query 的中间位置插入 150bp（测试大间隙处理）
        query.insert(500, generateRandomDNA(150, 5001));

        // 创建锚点：覆盖插入前后的区域
        anchor::Anchors anchors;

        // 锚点 1: 前半部分 (ref: 0-99, query: 0-99)
        anchor::Anchor a1;
        a1.hash = 3001;
        a1.rid_ref = 0;
        a1.pos_ref = 0;
        a1.rid_qry = 0;
        a1.pos_qry = 0;
        a1.span = 100;
        a1.is_rev = false;
        anchors.push_back(a1);

        // 锚点 2: 后半部分 (ref: 500-599, query: 650-749)
        // query 位置偏移了 150bp（插入长度）
        anchor::Anchor a2;
        a2.hash = 3002;
        a2.rid_ref = 0;
        a2.pos_ref = 500;
        a2.rid_qry = 0;
        a2.pos_qry = 650;
        a2.span = 100;
        a2.is_rev = false;
        anchors.push_back(a2);

        auto cigar = align::globalAlignMM2(ref, query, anchors);
        REQUIRE(cigar.size() > 0);

        std::string cigar_str = cigarToString(cigar);
        MESSAGE("CIGAR (large gap): ", cigar_str);

        // 验证完整覆盖
        std::size_t ref_len = cigar::getRefLength(cigar);
        std::size_t qry_len = cigar::getQueryLength(cigar);
        CHECK(ref_len == ref.size());
        CHECK(qry_len == query.size());
    }

    TEST_CASE("globalAlignMM2 - 与 globalAlignKSW2 结果一致性（无锚点）") {
        std::string ref = generateRandomDNA(500, 6000);
        std::string query = mutateSequence(ref, 0.02, 0.01, 6001);

        anchor::Anchors empty_anchors;
        auto cigar_mm2 = align::globalAlignMM2(ref, query, empty_anchors);
        auto cigar_ksw2 = align::globalAlignKSW2(ref, query);

        // 两者应该产生相同的结果（或至少长度相同）
        std::size_t mm2_ref_len = cigar::getRefLength(cigar_mm2);
        std::size_t mm2_qry_len = cigar::getQueryLength(cigar_mm2);
        std::size_t ksw2_ref_len = cigar::getRefLength(cigar_ksw2);
        std::size_t ksw2_qry_len = cigar::getQueryLength(cigar_ksw2);

        CHECK(mm2_ref_len == ksw2_ref_len);
        CHECK(mm2_qry_len == ksw2_qry_len);
        CHECK(mm2_ref_len == ref.size());
        CHECK(mm2_qry_len == query.size());
    }

    TEST_CASE("removeRefGapColumns - drop ref gap columns from aligned seq") {
        // 说明：测试"按 ref_gap_pos 删列"的纯过滤功能（原地修改）。
        // 输入序列已经是对齐后的（含 gap），本函数只负责删除参考为 gap 的列。

        // 输入：已对齐序列 "AC-GT"（长度 5）
        std::string seq = "AC-GT";

        // ref_gap_pos: true 表示参考该列是 gap，需要删除该列。
        // 我们删除第2列(0-based==2)这一列，输出应为 "ACGT"
        const std::vector<bool> ref_gap_pos = {false, false, true, false, false};

        align::RefAligner::removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "ACGT");
    }

    TEST_CASE("removeRefGapColumns - keep existing '-' as base when not in ref gap pos") {
        // 说明：输入序列本身含有 '-'，但只要 ref_gap_pos 该列为 false，就保留。

        std::string seq = "A-CG"; // 包含 1 个原生 '-'
        const std::vector<bool> ref_gap_pos = {false, false, false, false}; // 所有列都保留

        align::RefAligner::removeRefGapColumns(seq, ref_gap_pos);
        CHECK(seq == "A-CG");
    }


}

// ==================================================================
// 比对准确性测试套件
// ==================================================================
TEST_SUITE("align") {

    // ------------------------------------------------------------------
    // 辅助函数：验证 CIGAR 的正确性
    // ------------------------------------------------------------------
    static bool verifyCigar(const std::string& ref, const std::string& query,
                           const cigar::Cigar_t& cigar) {
        std::size_t ref_pos = 0;
        std::size_t qry_pos = 0;

        for (const auto& op : cigar) {
            char op_char;
            uint32_t len;
            cigar::intToCigar(op, op_char, len);

            switch (op_char) {
                case 'M':
                case '=':
                case 'X':
                    ref_pos += len;
                    qry_pos += len;
                    break;
                case 'I':
                    qry_pos += len;
                    break;
                case 'D':
                case 'N':
                    ref_pos += len;
                    break;
                case 'S':
                case 'H':
                    qry_pos += len;
                    break;
                default:
                    return false;
            }
        }

        return ref_pos == ref.size() && qry_pos == query.size();
    }

    // ------------------------------------------------------------------
    // 辅助函数：计算 CIGAR 的编辑距离（简化版）
    // ------------------------------------------------------------------
    static size_t getCigarEditDistance(const cigar::Cigar_t& cigar) {
        size_t edit_dist = 0;

        for (const auto& op : cigar) {
            char op_char;
            uint32_t len;
            cigar::intToCigar(op, op_char, len);

            if (op_char == 'I' || op_char == 'D') {
                edit_dist += len;
            }
        }

        return edit_dist;
    }

    // ------------------------------------------------------------------
    // 测试：高相似度（95%-99%）比对准确性
    // ------------------------------------------------------------------
    TEST_CASE("Accuracy - High similarity (95%-99%)") {
        constexpr int NUM_TESTS = 50;
        constexpr size_t SEQ_LEN = 5000;

        std::cout << "\n========== 高相似度比对准确性测试 (500bp, " << NUM_TESTS << " 次) ==========\n";

        struct AccuracyStats {
            size_t total_tests = 0;
            size_t valid_cigars = 0;
            size_t perfect_match = 0;
            double avg_edit_dist = 0.0;
        };

        std::vector<double> similarity_levels = {0.99, 0.98, 0.97, 0.95};

        for (double similarity : similarity_levels) {
            double snp_rate = (1.0 - similarity) * 0.7;
            double indel_rate = (1.0 - similarity) * 0.3;

            std::map<std::string, AccuracyStats> algorithm_stats;
            algorithm_stats["KSW2"];
            algorithm_stats["WFA2"];
            algorithm_stats["MM2"];

            // 锚点统计
            size_t total_anchors = 0;
            size_t min_anchors = std::numeric_limits<size_t>::max();
            size_t max_anchors = 0;

            for (int i = 0; i < NUM_TESTS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
                std::string query = mutateSequence(ref, snp_rate, indel_rate, i * 3 + 1);

                // 生成真实的锚点（基于 minimizer 匹配）
                // k=15, w=10：适合高相似度场景的参数
                anchor::Anchors anchors = generateRealAnchors(ref, query, 15, 10);

                // 统计锚点
                total_anchors += anchors.size();
                min_anchors = std::min(min_anchors, anchors.size());
                max_anchors = std::max(max_anchors, anchors.size());

                // KSW2
                {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    auto& stats = algorithm_stats["KSW2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }

                // WFA2
                {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    auto& stats = algorithm_stats["WFA2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }

                // MM2
                {
                    auto cigar = align::globalAlignMM2(ref, query, anchors);
                    auto& stats = algorithm_stats["MM2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        size_t edit_dist = getCigarEditDistance(cigar);
                        stats.avg_edit_dist += edit_dist;
                        if (edit_dist == 0) stats.perfect_match++;
                    }
                }
            }

            std::cout << "\n相似度 " << (similarity * 100) << "%:\n";
            std::cout << "  [锚点] 平均=" << (total_anchors / NUM_TESTS)
                      << ", 最小=" << (min_anchors == std::numeric_limits<size_t>::max() ? 0 : min_anchors)
                      << ", 最大=" << max_anchors << "\n";
            for (const auto& [name, stats] : algorithm_stats) {
                double accuracy = 100.0 * stats.valid_cigars / stats.total_tests;
                double avg_dist = stats.avg_edit_dist / stats.total_tests;
                std::cout << "  " << std::setw(6) << name
                          << ": 准确率=" << std::fixed << std::setprecision(1) << accuracy << "%"
                          << ", 平均编辑距离=" << std::setprecision(2) << avg_dist
                          << ", 完美匹配=" << stats.perfect_match << "/" << stats.total_tests << "\n";
            }
        }

        std::cout << "========================================================\n\n";
    }

    // ------------------------------------------------------------------
    // 测试：低相似度（70%-90%）比对准确性
    // ------------------------------------------------------------------
    TEST_CASE("Accuracy - Low similarity (70%-90%)") {
        constexpr int NUM_TESTS = 50;
        constexpr size_t SEQ_LEN = 5000;

        std::cout << "\n========== 低相似度比对准确性测试 (500bp, " << NUM_TESTS << " 次) ==========\n";

        struct AccuracyStats {
            size_t total_tests = 0;
            size_t valid_cigars = 0;
            size_t failed = 0;
            double avg_edit_dist = 0.0;
        };

        std::vector<double> similarity_levels = {0.90, 0.85, 0.80, 0.75, 0.70};

        for (double similarity : similarity_levels) {
            double snp_rate = (1.0 - similarity) * 0.6;
            double indel_rate = (1.0 - similarity) * 0.4;

            std::map<std::string, AccuracyStats> algorithm_stats;
            algorithm_stats["KSW2"];
            algorithm_stats["WFA2"];
            algorithm_stats["MM2"];

            for (int i = 0; i < NUM_TESTS; ++i) {
                std::string ref = generateRandomDNA(SEQ_LEN, i * 3);
                std::string query = mutateSequence(ref, snp_rate, indel_rate, i * 3 + 1);

                // 生成真实的锚点（低相似度场景：使用更小的 k 和 w）
                // k=13, w=8：低相似度时需要更宽松的参数以获得足够的锚点
                anchor::Anchors anchors = generateRealAnchors(ref, query, 13, 8);

                // KSW2
                {
                    auto cigar = align::globalAlignKSW2(ref, query);
                    auto& stats = algorithm_stats["KSW2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }

                // WFA2
                {
                    auto cigar = align::globalAlignWFA2(ref, query);
                    auto& stats = algorithm_stats["WFA2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }

                // MM2
                {
                    auto cigar = align::globalAlignMM2(ref, query, anchors);
                    auto& stats = algorithm_stats["MM2"];
                    stats.total_tests++;
                    if (verifyCigar(ref, query, cigar)) {
                        stats.valid_cigars++;
                        stats.avg_edit_dist += getCigarEditDistance(cigar);
                    } else {
                        stats.failed++;
                    }
                }
            }

            std::cout << "\n相似度 " << (similarity * 100) << "%:\n";
            for (const auto& [name, stats] : algorithm_stats) {
                double accuracy = 100.0 * stats.valid_cigars / stats.total_tests;
                double avg_dist = stats.avg_edit_dist / (stats.valid_cigars > 0 ? stats.valid_cigars : 1);
                std::cout << "  " << std::setw(6) << name
                          << ": 准确率=" << std::fixed << std::setprecision(1) << accuracy << "%"
                          << ", 平均编辑距离=" << std::setprecision(2) << avg_dist
                          << ", 失败次数=" << stats.failed << "/" << stats.total_tests << "\n";
            }
        }

        std::cout << "\n备注：低相似度时，某些算法可能因为参数限制而失败\n";
        std::cout << "========================================================\n\n";
    }



}

// ------------------------------------------------------------------
// 一键运行性能测试的辅助说明
// ------------------------------------------------------------------
// 说明：
// 1. 正确性测试默认运行，用于验证基本功能
// 2. 性能测试默认跳过（doctest::skip(true)），需要手动启用
// 3. 运行性能测试：
//    ./halign4_tests -tc="*Performance*" --no-skip
// 4. 只运行特定性能测试：
//    ./halign4_tests -tc="*Short*" --no-skip
// ------------------------------------------------------------------
