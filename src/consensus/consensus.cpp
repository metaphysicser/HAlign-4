#include "consensus.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <cereal/archives/json.hpp>

#if __has_include(<omp.h>)
    #include <omp.h>
#endif

namespace consensus
{
    // 选择共识碱基：仅在 A/C/G/T/U 之间选择（不会返回 N 或 '-'）。
    // 说明：如果这五者计数均为 0（例如该位点在输入中全为 N 或 gap），
    // 则按优先级回退为 'A'（可根据需要调整为抛出异常或其它策略）。
    char pickConsensusChar(const SiteCount& sc)
    {
        // 固定优先级（在相等时使用）：A > C > G > T > U
        std::uint32_t best = sc.a;
        char best_ch = 'A';

        auto upd = [&](std::uint32_t v, char ch) {
            if (v > best) { best = v; best_ch = ch; }
        };

        // 只比较 A/C/G/T/U，避免选出 N 或 gap
        upd(sc.c, 'C');
        upd(sc.g, 'G');
        upd(sc.t, 'T');
        upd(sc.u, 'U');

        return best_ch;
    }

    char pickConsensusCharWithGap(const SiteCount& sc)
    {
        const std::uint32_t best_non_gap = std::max({
            sc.a, sc.c, sc.g, sc.t, sc.u, sc.n
        });
        if (sc.dash > best_non_gap) {
            return '-';
        }
        return pickConsensusChar(sc);
    }

    // 将 consensus 序列按 FASTA 格式写出，行宽固定为 80
    // 注意：写入前会确保父目录存在
    void writeConsensusFasta(const FilePath& out_fasta, const std::string& seq)
    {
        file_io::ensureParentDirExists(out_fasta);

        std::ofstream ofs(out_fasta, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open fasta output: " + out_fasta.string());
        }

        ofs << ">consensus\n";
        constexpr std::size_t width = 80;
        for (std::size_t i = 0; i < seq.size(); i += width) {
            const std::size_t n = std::min(width, seq.size() - i);
            ofs.write(seq.data() + (std::streamoff)i, (std::streamsize)n);
            ofs.put('\n');
        }
    }

    // 使用 cereal 将计数写为 JSON（项目中 prefer cereal）
    void writeCountsJson(const FilePath& out_json, const ConsensusJson& cj)
    {
        file_io::ensureParentDirExists(out_json);

        std::ofstream ofs(out_json, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("failed to open json output: " + out_json.string());
        }

        cereal::JSONOutputArchive ar(ofs);
        ar(cereal::make_nvp("consensus", cj));
    }



    // 将单条序列的按列计数逻辑提取为一个独立函数。
    // 该函数只负责把单条序列的每个位点计数累加到 cj.counts 中，
    // 不负责更新 cj.num_seqs（由调用者在安全的语义下进行递增）。
    //
    // 并行化策略：对单条序列，在 "位置" 维度并行（每个线程负责不同的列索引），
    // 因为每个列对应的 SiteCount 是独立的，多个线程不会写入同一索引，故无需原子操作。
    // 注意：如果并行处理多条序列（在调用层面并行），必须保证不同线程不会并发写入同一列。
    //
    // 详细优化说明（中文，供维护者阅读）
    //
    // 总体目标：
    // - 在处理大量对齐序列时，将“按列统计碱基计数”这个瓶颈路径尽量并行化和向量化，
    //   同时避免高成本的原子操作或频繁的锁竞争，尽可能利用缓存局部性与 SIMD 指令。
    //
    // 主要优化手段（本文件实现中可见）：
    // 1) 分支最小化（branch minimization）
    //    - 传统的 switch/case 每个字符要走分支，分支预测失败开销大。
    //    - 用 (idx == constant) 转成 0/1 的加法消除分支，编译器更容易生成向量指令。
    //
    // 2) 循环展开（loop unrolling）与预取（prefetch）
    //    - 通过 4-路（或可调整）循环展开减少循环控制开销、增加指令级并行（ILP）。
    //    - 使用 __builtin_prefetch 提前载入将要写入的缓存行，降低缓存未命中延迟（对长序列/大工作集有效）。
    //
    // 3) 线程局部（thread-local）累加 + 批处理（batching）
    //    - 直接并发写全局 counts 会导致 cache-line 冲突（false sharing）。为避免，采用每个线程维护本地计数数组，
    //      在处理完一批序列后再按列合并（reduce）。这种策略通过牺牲一些内存来大幅降低并发写冲突，
    //      在线程数与对齐长度适中时通常能得到最好吞吐。
    //    - 批大小（batch_size）需要调优：过小会导致调度/同步开销，过大会占用更多内存且增加延迟。
    //
    // 4) SoA (Structure of Arrays) vs AoS (Array of Structures)
    //    - AoS（SiteCount 存在 cj.counts[i]）在合并或向量化时不便，因为每个 SiteCount 中字段交错。
    //    - SoA 将每个碱基的计数放入独立的数组（A[], C[], G[]...），便于对单个碱基列进行连续读写，
    //      更有利于 SIMD 与缓存预取。实现中提供了 SoA 的批处理路径以求最大化性能。
    //
    // 5) OpenMP + simd 指示（#pragma omp parallel / #pragma omp simd / ivdep）
    //    - 使用 OpenMP 做线程级并行（按位置或按序列分配工作）；同时在内循环增加 simd 提示以帮助编译器向量化。
    //
    // 6) 编译器与编译选项
    //    - 在 CMake 中对 Release 模式启用 -O3、-march=native、可选 -flto，有助于生成高效的向量指令与内联。
    //
    // 适用场景与权衡：
    // - 当 aln_len（对齐长度）很大（几千到几万）且每批序列数量也大时，SoA + 线程本地合并通常最快。
    // - 当 aln_len 很小（例如 < 256）时，线程并行化的开销可能抵消收益，应退化到单线程或较小线程数。
    // - 本地计数会使用额外内存：大线程数与大对齐长度会带来显著内存占用（T * aln_len * sizeof(count)）。
    // - 合并阶段仍然要按列扫描，一次性合并成本需在批大小与内存之间权衡。
    //
    // 实践建议：
    // - 在目标机器上做小规模基准（不同 batch_size、线程数、aln_len）来选择最优参数。
    // - 结合性能分析工具（perf / VTune / likwid）观察缓存未命中、内存带宽与分支失误。
    // - 若需要极限性能，可进一步用 AVX2/AVX512 intrinsics 对 SoA 路径做微调。
    static void processSequenceParallel(const std::string& s, ConsensusJson& cj, int thread)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (s.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: expect " + std::to_string(aln_len) +
                                     ", got " + std::to_string(s.size()));
        }

        // 预取常量以减少循环内查找开销
        const unsigned char* data = reinterpret_cast<const unsigned char*>(s.data());
        std::uint8_t* base_map = const_cast<std::uint8_t*>(k_base_map.data());
        SiteCount* counts = cj.counts.data();

        const std::size_t limit = (aln_len / 4) * 4;

#if __has_include(<omp.h>)
        // 使用 parallel for simd，允许 OpenMP 分配线程并启用向量化
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = 0; i < limit; i += 4) {
            // 轻量预取未来缓存行，距离可调（16~64 bytes => 16 positions 粗略）
            __builtin_prefetch(&counts[i + 16]);

            const std::uint8_t idx0 = base_map[data[i + 0]];
            const std::uint8_t idx1 = base_map[data[i + 1]];
            const std::uint8_t idx2 = base_map[data[i + 2]];
            const std::uint8_t idx3 = base_map[data[i + 3]];

            SiteCount& sc0 = counts[i + 0];
            SiteCount& sc1 = counts[i + 1];
            SiteCount& sc2 = counts[i + 2];
            SiteCount& sc3 = counts[i + 3];

            sc0.a += static_cast<std::uint32_t>(idx0 == 0);
            sc0.c += static_cast<std::uint32_t>(idx0 == 1);
            sc0.g += static_cast<std::uint32_t>(idx0 == 2);
            sc0.t += static_cast<std::uint32_t>(idx0 == 3);
            sc0.u += static_cast<std::uint32_t>(idx0 == 4);
            sc0.n += static_cast<std::uint32_t>(idx0 == 5);
            sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

            sc1.a += static_cast<std::uint32_t>(idx1 == 0);
            sc1.c += static_cast<std::uint32_t>(idx1 == 1);
            sc1.g += static_cast<std::uint32_t>(idx1 == 2);
            sc1.t += static_cast<std::uint32_t>(idx1 == 3);
            sc1.u += static_cast<std::uint32_t>(idx1 == 4);
            sc1.n += static_cast<std::uint32_t>(idx1 == 5);
            sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

            sc2.a += static_cast<std::uint32_t>(idx2 == 0);
            sc2.c += static_cast<std::uint32_t>(idx2 == 1);
            sc2.g += static_cast<std::uint32_t>(idx2 == 2);
            sc2.t += static_cast<std::uint32_t>(idx2 == 3);
            sc2.u += static_cast<std::uint32_t>(idx2 == 4);
            sc2.n += static_cast<std::uint32_t>(idx2 == 5);
            sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

            sc3.a += static_cast<std::uint32_t>(idx3 == 0);
            sc3.c += static_cast<std::uint32_t>(idx3 == 1);
            sc3.g += static_cast<std::uint32_t>(idx3 == 2);
            sc3.t += static_cast<std::uint32_t>(idx3 == 3);
            sc3.u += static_cast<std::uint32_t>(idx3 == 4);
            sc3.n += static_cast<std::uint32_t>(idx3 == 5);
            sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
        }

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(thread)
#endif
        for (std::size_t i = limit; i < aln_len; ++i) {
            const std::uint8_t idx = base_map[data[i]];
            SiteCount& sc = counts[i];
            sc.a += static_cast<std::uint32_t>(idx == 0);
            sc.c += static_cast<std::uint32_t>(idx == 1);
            sc.g += static_cast<std::uint32_t>(idx == 2);
            sc.t += static_cast<std::uint32_t>(idx == 3);
            sc.u += static_cast<std::uint32_t>(idx == 4);
            sc.n += static_cast<std::uint32_t>(idx == 5);
            sc.dash += static_cast<std::uint32_t>(idx == 6);
        }
    }

    // 批量并行处理：每个线程维护自己的本地 counts 数组来累加多条序列，
    // 处理完批后再把本地 counts 合并到全局 cj.counts。这能显著减少对全局 counts 的并发写入，
    // 降低 false sharing 并提高缓存局部性。适合 aln_len 和 batch_size 都较大的场景。
    static void processBatchParallel(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        // 先在主线程中验证所有序列长度一致，避免并行区抛异常不安全
        for (const auto& s : seqs) {
            if (s.size() != aln_len) {
                throw std::runtime_error("alignment length mismatch in batch processing");
            }
        }

        // 扁平化本地计数：一块连续内存，大小为 T * aln_len
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals.data() + (std::size_t)tid * aln_len;

            // 对单条序列按位置累加到线程本地计数；采用分支最小化的布尔加法
            for (std::size_t i = 0; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并本地 counts 到全局 cj.counts（按列并行）
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            // 使用 64-bit 临时累加以避免溢出（尽管 SiteCount 是 uint32）
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     下面对批处理（线程本地累加 + 合并）函数做详细注释：

     processBatchParallelWithLocals:
     - 每线程拥有一块连续的 SiteCount 数组（locals），大小为 aln_len。
     - 每个线程将它负责的序列累加到本地 locals（避免并发写入全局 cj.counts）。
     - 合并阶段按列并行：每个线程负责合并若干列，将本地累加的结果加到全局 counts。

     优点：
     - 极大降低 false sharing（不同线程不会频繁写同一 cache line）。
     - 简单实现，易于理解与调试。

     局限：
     - 每线程的 locals 使用较多内存（T * aln_len * sizeof(SiteCount)）。
     - SiteCount 内部字段为 AoS，合并时访问字段交错，向量化受限。
    */

    // 批处理（使用外部分配的 locals 缓冲以避免每批分配开销）
    static void processBatchParallelWithLocals(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads, std::vector<SiteCount>& locals)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();
        SiteCount* locals_ptr = locals.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());
            SiteCount* local = locals_ptr + (std::size_t)tid * aln_len;

            // 更高效的按位置更新：4-路展开并使用布尔加法，减少分支
            std::size_t i = 0;
            const std::size_t limit = (aln_len / 4) * 4;
            for (; i < limit; i += 4) {
                __builtin_prefetch(&local[i + 16]);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                SiteCount& sc0 = local[i + 0];
                SiteCount& sc1 = local[i + 1];
                SiteCount& sc2 = local[i + 2];
                SiteCount& sc3 = local[i + 3];

                sc0.a += static_cast<std::uint32_t>(idx0 == 0);
                sc0.c += static_cast<std::uint32_t>(idx0 == 1);
                sc0.g += static_cast<std::uint32_t>(idx0 == 2);
                sc0.t += static_cast<std::uint32_t>(idx0 == 3);
                sc0.u += static_cast<std::uint32_t>(idx0 == 4);
                sc0.n += static_cast<std::uint32_t>(idx0 == 5);
                sc0.dash += static_cast<std::uint32_t>(idx0 == 6);

                sc1.a += static_cast<std::uint32_t>(idx1 == 0);
                sc1.c += static_cast<std::uint32_t>(idx1 == 1);
                sc1.g += static_cast<std::uint32_t>(idx1 == 2);
                sc1.t += static_cast<std::uint32_t>(idx1 == 3);
                sc1.u += static_cast<std::uint32_t>(idx1 == 4);
                sc1.n += static_cast<std::uint32_t>(idx1 == 5);
                sc1.dash += static_cast<std::uint32_t>(idx1 == 6);

                sc2.a += static_cast<std::uint32_t>(idx2 == 0);
                sc2.c += static_cast<std::uint32_t>(idx2 == 1);
                sc2.g += static_cast<std::uint32_t>(idx2 == 2);
                sc2.t += static_cast<std::uint32_t>(idx2 == 3);
                sc2.u += static_cast<std::uint32_t>(idx2 == 4);
                sc2.n += static_cast<std::uint32_t>(idx2 == 5);
                sc2.dash += static_cast<std::uint32_t>(idx2 == 6);

                sc3.a += static_cast<std::uint32_t>(idx3 == 0);
                sc3.c += static_cast<std::uint32_t>(idx3 == 1);
                sc3.g += static_cast<std::uint32_t>(idx3 == 2);
                sc3.t += static_cast<std::uint32_t>(idx3 == 3);
                sc3.u += static_cast<std::uint32_t>(idx3 == 4);
                sc3.n += static_cast<std::uint32_t>(idx3 == 5);
                sc3.dash += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                SiteCount& sc = local[i];
                sc.a += static_cast<std::uint32_t>(idx == 0);
                sc.c += static_cast<std::uint32_t>(idx == 1);
                sc.g += static_cast<std::uint32_t>(idx == 2);
                sc.t += static_cast<std::uint32_t>(idx == 3);
                sc.u += static_cast<std::uint32_t>(idx == 4);
                sc.n += static_cast<std::uint32_t>(idx == 5);
                sc.dash += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并本地 counts 到全局 cj.counts（按列并行）
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            for (int t = 0; t < T; ++t) {
                const SiteCount& ls = locals_ptr[(std::size_t)t * aln_len + i];
                sa += ls.a; sc_ += ls.c; sg += ls.g; st += ls.t; su += ls.u; sn += ls.n; sd += ls.dash;
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }

    /*
     processBatchParallelWithSoA（SoA 版本）：
     - 对每个碱基维护独立数组（A[], C[], G[], T[], U[], N[], Dash[]），并为每个线程分配一段连续内存（T * aln_len）。
     - 在累加阶段，每线程对自己的段按位置累加：a_ptr[i]++ 等。由于同一种碱基的计数是连续的，
       读取/写入模式更友好于 CPU 的向量化与缓存预取策略。
     - 合并阶段按列（i）对每个线程的对应位置求和并写回 cj.counts[i]。
     *
     * 优势（为什么更快）：
     * 1) 向量化友好：SoA 允许编译器把一条指令应用到连续的 a_ptr[] 元素，生成 SIMD 指令（如 AVX2 之类）。
     * 2) 减少内存带宽浪费：当处理某个碱基时只触碰该碱基的数组，减少对其它字段不必要的缓存写回。
     * 3) 合并步骤中按列累加多个线程的连续内存，内存访问模式简单，有利于预取与硬件合并。
     *
     * 代价：
     * - 额外的内存开销（7 个 uint32_t 数组 * T * aln_len），但通常比原子/锁/频繁缓存同步更划算。
     * - 实现复杂度略增，但对性能敏感场景值得。
    */

    static void processBatchParallelWithSoA(const std::vector<std::string>& seqs, ConsensusJson& cj, int threads,
                                            std::vector<std::uint32_t>& localsA,
                                            std::vector<std::uint32_t>& localsC,
                                            std::vector<std::uint32_t>& localsG,
                                            std::vector<std::uint32_t>& localsT,
                                            std::vector<std::uint32_t>& localsU,
                                            std::vector<std::uint32_t>& localsN,
                                            std::vector<std::uint32_t>& localsDash)
    {
        const std::size_t aln_len = (std::size_t)cj.aln_len;
        if (aln_len == 0 || seqs.empty()) return;

        int T = threads;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        T = 1;
#endif

        const std::uint8_t* base_map = k_base_map.data();

#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = omp_get_thread_num();
#else
        for (std::size_t s = 0; s < seqs.size(); ++s) {
            const int tid = 0;
#endif
            const std::string& str = seqs[s];
            const unsigned char* data = reinterpret_cast<const unsigned char*>(str.data());

            const std::size_t base_off = (std::size_t)tid * aln_len;
            std::uint32_t* a_ptr = localsA.data() + base_off;
            std::uint32_t* c_ptr = localsC.data() + base_off;
            std::uint32_t* g_ptr = localsG.data() + base_off;
            std::uint32_t* t_ptr = localsT.data() + base_off;
            std::uint32_t* u_ptr = localsU.data() + base_off;
            std::uint32_t* n_ptr = localsN.data() + base_off;
            std::uint32_t* d_ptr = localsDash.data() + base_off;

            const std::size_t limit = (aln_len / 4) * 4;
            std::size_t i = 0;
            for (; i < limit; i += 4) {
                __builtin_prefetch(a_ptr + i + 16);
                const std::uint8_t idx0 = base_map[data[i + 0]];
                const std::uint8_t idx1 = base_map[data[i + 1]];
                const std::uint8_t idx2 = base_map[data[i + 2]];
                const std::uint8_t idx3 = base_map[data[i + 3]];

                a_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 0);
                c_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 1);
                g_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 2);
                t_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 3);
                u_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 4);
                n_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 5);
                d_ptr[i + 0] += static_cast<std::uint32_t>(idx0 == 6);

                a_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 0);
                c_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 1);
                g_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 2);
                t_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 3);
                u_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 4);
                n_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 5);
                d_ptr[i + 1] += static_cast<std::uint32_t>(idx1 == 6);

                a_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 0);
                c_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 1);
                g_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 2);
                t_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 3);
                u_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 4);
                n_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 5);
                d_ptr[i + 2] += static_cast<std::uint32_t>(idx2 == 6);

                a_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 0);
                c_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 1);
                g_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 2);
                t_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 3);
                u_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 4);
                n_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 5);
                d_ptr[i + 3] += static_cast<std::uint32_t>(idx3 == 6);
            }
            for (; i < aln_len; ++i) {
                const std::uint8_t idx = base_map[data[i]];
                a_ptr[i] += static_cast<std::uint32_t>(idx == 0);
                c_ptr[i] += static_cast<std::uint32_t>(idx == 1);
                g_ptr[i] += static_cast<std::uint32_t>(idx == 2);
                t_ptr[i] += static_cast<std::uint32_t>(idx == 3);
                u_ptr[i] += static_cast<std::uint32_t>(idx == 4);
                n_ptr[i] += static_cast<std::uint32_t>(idx == 5);
                d_ptr[i] += static_cast<std::uint32_t>(idx == 6);
            }
        }

        // 合并：按列读取 locals arrays 并写入 cj.counts
#if __has_include(<omp.h>)
        #pragma omp parallel for schedule(static) num_threads(T)
#endif
        for (std::size_t i = 0; i < aln_len; ++i) {
            std::uint64_t sa = 0, sc_ = 0, sg = 0, st = 0, su = 0, sn = 0, sd = 0;
            const std::size_t stride = aln_len;
            for (int t = 0; t < T; ++t) {
                const std::size_t off = (std::size_t)t * stride + i;
                sa += localsA[off]; sc_ += localsC[off]; sg += localsG[off]; st += localsT[off]; su += localsU[off]; sn += localsN[off]; sd += localsDash[off];
            }
            SiteCount& dst = cj.counts[i];
            dst.a += static_cast<std::uint32_t>(sa);
            dst.c += static_cast<std::uint32_t>(sc_);
            dst.g += static_cast<std::uint32_t>(sg);
            dst.t += static_cast<std::uint32_t>(st);
            dst.u += static_cast<std::uint32_t>(su);
            dst.n += static_cast<std::uint32_t>(sn);
            dst.dash += static_cast<std::uint32_t>(sd);
        }
    }


    ConsensusResult generateConsensusResult(const FilePath& aligned_fasta,
                                            const FilePath& out_fasta,
                                            const FilePath& out_json,
                                            std::uint64_t seq_limit,
                                            int thread,
                                            size_t batch_size)
    {
        file_io::requireRegularFile(aligned_fasta, "aligned_fasta");
        if (batch_size == 0) {
            batch_size = 4096;
        }

        seq_io::KseqReader reader(aligned_fasta);

        // 读第一条确定 aln_len
        seq_io::SeqRecord rec;
        if (!reader.next(rec)) {
            throw std::runtime_error("aligned fasta is empty: " + aligned_fasta.string());
        }

        const std::size_t aln_len = rec.seq.size();
        if (aln_len == 0) {
            throw std::runtime_error("first sequence length is 0: " + aligned_fasta.string());
        }

        ConsensusResult result;
        ConsensusJson& cj = result.counts;
        cj.aln_len = (std::uint64_t)aln_len;
        cj.counts.assign(aln_len, SiteCount{});

        std::uint64_t num_seqs = 0;

        // 批处理参数：批大小可调（经验值），在内存允许的范围内放大能减少调度开销
        std::vector<std::string> batch;
        batch.reserve(batch_size + 1);

        // 首条先放入 batch（并检查长度）
        if (rec.seq.size() != aln_len) {
            throw std::runtime_error("alignment length mismatch: first record length changed");
        }
        batch.push_back(std::move(rec.seq));
        ++num_seqs;

        // 预分配并复用 thread-local 缓冲
        int T = thread;
#if __has_include(<omp.h>)
        if (T <= 0) T = omp_get_max_threads();
        if (T <= 0) T = 1;
#else
        if (T <= 0) T = 1;
#endif
        std::vector<SiteCount> locals;
        try {
            locals.assign((std::size_t)T * aln_len, SiteCount{});
        } catch (...) {
            throw std::runtime_error("failed to allocate thread-local counts");
        }

        auto flush_batch = [&]() {
            if (batch.empty()) {
                return;
            }
            std::memset(locals.data(), 0, locals.size() * sizeof(SiteCount));
            processBatchParallelWithLocals(batch, cj, thread, locals);
            batch.clear();
        };

        // 读取并按批处理；注意不要在批次之间额外读取一条记录，避免跳过序列。
        while (seq_limit == 0 || num_seqs < seq_limit) {
            if (batch.size() >= batch_size) {
                flush_batch();
            }

            if (!reader.next(rec)) {
                break;
            }
            if (rec.seq.size() != aln_len) {
                throw std::runtime_error("alignment length mismatch when reading");
            }
            batch.push_back(std::move(rec.seq));
            ++num_seqs;
        }
        flush_batch();

        cj.num_seqs = num_seqs;

        if (num_seqs == 0) {
            throw std::runtime_error("no sequences processed");
        }

        // 生成共识序列（单线程选多数）
        result.gap_seq.reserve(aln_len);
        result.seq.reserve(aln_len);
        for (std::size_t i = 0; i < aln_len; ++i) {
            const char ch = pickConsensusCharWithGap(cj.counts[i]);
            result.gap_seq.push_back(ch);
            if (ch != '-') {
                result.seq.push_back(ch);
            }
        }

        writeConsensusFasta(out_fasta, result.seq);
        writeCountsJson(out_json, cj);

        return result;
    }

    std::string generateConsensusSequence(const FilePath& aligned_fasta,
                                                       const FilePath& out_fasta,
                                                       const FilePath& out_json,
                                                       std::uint64_t seq_limit,
                                                       int thread,
                                                       size_t batch_size)
    {
        return generateConsensusResult(
            aligned_fasta, out_fasta, out_json, seq_limit, thread, batch_size).seq;
    }


} // namespace consensus
