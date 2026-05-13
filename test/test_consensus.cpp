#include <doctest/doctest.h>
#include <thread>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
    #include <sys/resource.h>
#endif

#include "consensus.h"

namespace fs = std::filesystem;


// ------------------------- correctness helpers -------------------------

static void writeTextFile(const fs::path& p, const std::string& content) {
    std::ofstream ofs(p, std::ios::binary);
    REQUIRE_MESSAGE(ofs.good(), "cannot write file: " << p.string());
    ofs << content;
}

static std::string readSingleFastaSequence(const fs::path& fasta) {
    std::ifstream ifs(fasta, std::ios::binary);
    REQUIRE_MESSAGE(ifs.good(), "cannot read fasta: " << fasta.string());

    std::string line, seq;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        seq += line;
    }
    return seq;
}

static fs::path makeTempDirSimple() {
    fs::path dir = fs::temp_directory_path() / "halign4_tests_simple";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    REQUIRE_MESSAGE(!ec, "cannot create temp dir: " << dir.string() << " (" << ec.message() << ")");
    return dir;
}

// ------------------------- perf helpers -------------------------

// ------------- env helpers -------------

static std::size_t envSizeT(const char* name, std::size_t def) {
    if (const char* v = std::getenv(name); v && *v) {
        try { return (std::size_t)std::stoull(v); } catch (...) {}
    }
    return def;
}

static int envInt(const char* name, int def) {
    if (const char* v = std::getenv(name); v && *v) {
        try { return std::stoi(v); } catch (...) {}
    }
    return def;
}

// perf 开关：你已有
static bool perfEnabled() {
    const char* v = std::getenv("HALIGN4_RUN_PERF");
    return (v != nullptr) && (*v != '\0') && (std::string(v) != "0");
}
static bool shouldSkipPerf() { return !perfEnabled(); }

// perf 最大 n：默认 100k，避免默认就尝试 1e6
static std::size_t perfMaxN() {
    return envSizeT("HALIGN4_PERF_MAX_N", 100'000);
}

// perf 目录：大数据强制要求用户显式指定，避免写到 build 目录或 /tmp
static fs::path perfBaseDir() {
    if (const char* p = std::getenv("HALIGN4_PERF_DIR"); p && *p) {
        return fs::path(p);
    }
    return fs::path(); // empty => 未指定
}

// 估算 FASTA 大小（非常粗略，但足够做预检）
static std::uint64_t estimateFastaBytes(std::size_t n_seqs, std::size_t len) {
    // header: ">s<idx>\n" 约 12~20 字节；我们按 24 估
    // seq: len + "\n"
    const std::uint64_t per = (std::uint64_t)len + 1ull + 24ull;
    return per * (std::uint64_t)n_seqs;
}

static bool hasEnoughDisk(const fs::path& dir, std::uint64_t need_bytes) {
    std::error_code ec;
    auto sp = fs::space(dir, ec);
    if (ec) return true; // 取不到就不拦（但会打印提示）
    // 留 20% buffer，避免写到一半触发配额/满盘
    return sp.available > (need_bytes * 12ull) / 10ull;
}

// ------------- perf io -------------

static void writeLargeAlignedFasta(const fs::path& p, std::size_t n_seqs, std::size_t len) {
    std::ofstream ofs(p, std::ios::binary);
    ofs.exceptions(std::ios::badbit | std::ios::failbit);
    REQUIRE_MESSAGE(ofs.good(), "cannot write fasta: " << p.string());

    std::string seq(len, 'A');
    seq.push_back('\n');

    // 可选：写入进度（避免外部 watchdog 认为“无输出卡死”）
    const std::size_t progress_step = envSizeT("HALIGN4_PERF_PROGRESS_STEP", 0);

    for (std::size_t i = 0; i < n_seqs; ++i) {
        char hdr[64];
        int hn = std::snprintf(hdr, sizeof(hdr), ">s%zu\n", i);
        ofs.write(hdr, hn);
        ofs.write(seq.data(), (std::streamsize)seq.size());

        if (progress_step && (i % progress_step == 0)) {
            MESSAGE("writing fasta... i=" << i << "/" << n_seqs);
        }
    }
}

// ------------- perf runner -------------

static void runOnePerf(std::size_t n_seqs) {
    constexpr std::size_t LEN   = 30000;

    // 简化版：仅在显式开启 perf 时运行；写入输入（如需），然后直接调用共识计算并计时。
    if (!perfEnabled()) {
        MESSAGE("perf disabled - skip");
        return;
    }

    fs::path base = perfBaseDir();
    if (base.empty()) base = fs::temp_directory_path();

    fs::path dir = base / "halign4_tests_perf";
    std::error_code ec;
    fs::create_directories(dir, ec);
    REQUIRE_MESSAGE(!ec, "cannot create perf dir: " << dir.string() << " (" << ec.message() << ")");

    fs::path in_fa  = dir / ("aligned_" + std::to_string(n_seqs) + ".fasta");
    fs::path out_fa = dir / ("cons_"    + std::to_string(n_seqs) + ".fasta");
    fs::path out_js = dir / ("counts_"  + std::to_string(n_seqs) + ".json");

    const bool reuse = envInt("HALIGN4_PERF_REUSE_INPUT", 1) != 0;
    int threads_override = envInt("HALIGN4_PERF_THREADS", 0);
    if (threads_override <= 0) {
        unsigned int hc = std::thread::hardware_concurrency();
        threads_override = static_cast<int>(hc ? hc : 1u);
    }

    auto cleanup = [&]() {
        std::error_code e2;
        fs::remove(in_fa, e2);
        fs::remove(out_fa, e2);
        fs::remove(out_js, e2);
    };

    try {
        MESSAGE("perf start: n_seqs=" << n_seqs << " len=" << LEN << " threads=" << threads_override);

        if (!reuse || !fs::exists(in_fa)) {
            writeLargeAlignedFasta(in_fa, n_seqs, LEN);
        } else {
            MESSAGE("reuse existing input: " << in_fa.string());
        }

        auto t0 = std::chrono::steady_clock::now();
        (void)consensus::generateConsensusSequence(
            in_fa, out_fa, out_js,
            /*seq_limit=*/0, threads_override
        );
        auto t1 = std::chrono::steady_clock::now();

        MESSAGE("compute_done: seconds=" << std::chrono::duration<double>(t1 - t0).count()
                << " out_fa=" << out_fa.string() << " out_js=" << out_js.string());

        if (envInt("HALIGN4_PERF_CLEANUP", 1) != 0) {
            cleanup();
        }
    } catch (...) {
        if (envInt("HALIGN4_PERF_CLEANUP", 1) != 0) {
            cleanup();
        }
        throw;
    }
}


// ========================= correctness suite =========================

TEST_SUITE("consensus") {

TEST_CASE("generateConsensusResult - keeps gap-majority consensus and writes ungapped consensus") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned.fasta";
    fs::path out_fa = dir / "consensus.fasta";
    fs::path out_js = dir / "counts.json";

    const std::string aligned =
        ">s1\nACGT-\n"
        ">s2\nAC-T-\n"
        ">s3\nACGT-\n";

    writeTextFile(in_fa, aligned);

    consensus::ConsensusResult result = consensus::generateConsensusResult(
        in_fa, out_fa, out_js,
        0, 1
    );

    CHECK(result.gap_seq == "ACGT-");
    CHECK(result.seq == "ACGT");
    CHECK(readSingleFastaSequence(out_fa) == "ACGT");

    std::error_code ec;
    CHECK(fs::exists(out_js, ec));
    if (!ec && fs::exists(out_js, ec)) {
        CHECK(fs::file_size(out_js, ec) > 0);
    }
}

TEST_CASE("generateConsensusSequence - tie breaks to A (A > C > G > T > U)") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned_tie.fasta";
    fs::path out_fa = dir / "consensus_tie.fasta";
    fs::path out_js = dir / "counts_tie.json";

    const std::string aligned =
        ">s1\nA\n"
        ">s2\nC\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1
    );

    CHECK(cons == "A");
    CHECK(readSingleFastaSequence(out_fa) == "A");
}

TEST_CASE("generateConsensusResult - all gaps stay only in gap consensus") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "aligned_allgap.fasta";
    fs::path out_fa = dir / "consensus_allgap.fasta";
    fs::path out_js = dir / "counts_allgap.json";

    const std::string aligned =
        ">s1\n---\n"
        ">s2\n---\n"
        ">s3\n---\n";

    writeTextFile(in_fa, aligned);

    consensus::ConsensusResult result = consensus::generateConsensusResult(
        in_fa, out_fa, out_js,
        0, 1
    );

    CHECK(result.gap_seq == "---");
    CHECK(result.seq.empty());
    CHECK(readSingleFastaSequence(out_fa).empty());
}

TEST_CASE("generateConsensusSequence - single sequence returns same") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "single.fasta";
    fs::path out_fa = dir / "single_cons.fasta";
    fs::path out_js = dir / "single_counts.json";

    const std::string seq = ">s1\nACGTACGT\n";
    writeTextFile(in_fa, seq);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1
    );

    CHECK(cons == "ACGTACGT");
    CHECK(readSingleFastaSequence(out_fa) == "ACGTACGT");
}

TEST_CASE("generateConsensusSequence - U (uracil) handling") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "u_vs_t.fasta";
    fs::path out_fa = dir / "u_vs_t_cons.fasta";
    fs::path out_js = dir / "u_vs_t_counts.json";

    const std::string aligned =
        ">s1\nUAAAA\n"
        ">s2\nTAAAA\n"
        ">s3\nUAAAA\n";

    writeTextFile(in_fa, aligned);

    std::string cons = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1
    );

    CHECK(cons.size() == 5);
    CHECK(cons[0] == 'U');
}

TEST_CASE("generateConsensusSequence - seq_limit affects result") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "limit.fasta";
    fs::path out_fa = dir / "limit_cons.fasta";
    fs::path out_js = dir / "limit_counts.json";

    std::string content;
    content += ">s0\nA\n";
    content += ">s1\nA\n";
    content += ">s2\nC\n";
    content += ">s3\nC\n";
    content += ">s4\nC\n";
    writeTextFile(in_fa, content);

    std::string cons_all = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        0, 1
    );
    CHECK(cons_all == "C");

    std::string cons_lim = consensus::generateConsensusSequence(
        in_fa, out_fa, out_js,
        2, 1
    );
    CHECK(cons_lim == "A");
}

TEST_CASE("generateConsensusSequence - empty input throws") {
    auto dir = makeTempDirSimple();
    fs::path in_fa  = dir / "empty.fasta";
    fs::path out_fa = dir / "empty_cons.fasta";
    fs::path out_js = dir / "empty_counts.json";

    writeTextFile(in_fa, "");

    CHECK_THROWS_AS(
        consensus::generateConsensusSequence(in_fa, out_fa, out_js, 0, 1),
        std::runtime_error
    );
}

} // TEST_SUITE(consensus_generate)

// ========================= perf suite (opt-in) =========================

TEST_SUITE("consensus" * doctest::skip(shouldSkipPerf())) {

TEST_CASE("len=30000 n=100 batch=512 threads=max") {
    runOnePerf(100);
}


    TEST_CASE("len=30000 n=1000 batch=512 threads=max") {
    runOnePerf(1000);
}

TEST_CASE("len=30000 n=10000 batch=512 threads=max") {
    runOnePerf(10000);
}
//
TEST_CASE("len=30000 n=100000 batch=512 threads=max") {
     runOnePerf(100000);
}

} // TEST_SUITE(consensus_perf)
