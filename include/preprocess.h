#ifndef HALIGN4_PREPROCESS_H
#define HALIGN4_PREPROCESS_H

#include <cstddef>
#include "config.hpp"
#include "utils.h"
#include "consensus.h"
#include <algorithm>
#include <cctype>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ==============================================================
// 预处理模块（preprocess）头文件说明（详细中文注释）
//
// 本模块负责将用户输入的原始 FASTA（可以是本地路径或远程 URL）准备为后续分析的标准化数据，
// 并为共识计算/比对准备必要的中间文件。主要职责包括：
//  1) 将输入文件复制或下载到工作目录下的 `data/raw`；
//  2) 逐条读取输入序列并做“清洗/规范化”（例如大写化、把非 A/C/G/T/U 替换为 N，去除非法字符等）；
//  3) 将清洗后的序列写入 `data/clean`；
//  4) 维护一个 Top-K 选择器（长度优先）来挑选用于构建共识的候选序列集合（写入 `consensus_unaligned.fasta`）；
//  5) 返回处理的序列总数，供上层决定是否需要后续合并/更多处理。
//
// 重要语义与约定：
// - `workdir`：调用方提供的工作目录路径；本模块会在该目录下创建必要的子目录（data/raw, data/clean 等），
//   若目录不存在将尝试创建；如果要求为空（由上层传入并检查），则会在空目录中创建数据结构；
// - I/O 行为：如果 `input_path` 是远程 URL（例如 http(s):// 或以 // 开头），本模块会下载到本地；
//   否则会拷贝本地文件到工作目录。下载/拷贝失败会抛出异常（std::runtime_error）。
// - 异常与错误处理：函数在遇到严重 I/O 或解析错误时会抛出异常（std::runtime_error）；上层应捕获并记录。
// - 返回值：函数返回处理的记录数量（uint_t），若数量超过项目配置上限（config.hpp 中定义的 U_MAX），
//   值会被截断为 U_MAX 并记录告警（调用者应注意）。
//
// 性能与并发注意事项：
// - 该函数为 I/O 密集型：对大文件（GB 级）应关注磁盘带宽与缓冲（可通过 utils::seq_io 的 io 缓冲调优）；
// - 在高并发环境下，不要并行调用本函数写入同一 `workdir`，以免出现竞态；若需并行，使用不同的工作目录或外部协调。
//
// ==============================================================

// 预处理输入 FASTA，并返回处理的序列数量（total records processed）
//
// 参数：
//  - input_path: 输入 FASTA 文件路径，支持本地路径或远程 URL（字符串）。
//  - workdir:   工作目录（字符串），会在该目录下创建 data/raw 和 data/clean 子目录并写入中间文件。
//  - cons_n:    要为后续共识选择的序列数量（Top-K，按序列长度挑选），默认为 1000。
//
// 返回值：
//  - 返回实际处理的序列数（uint_t）；若处理条目过多超过 U_MAX，会截断为 U_MAX 并记录警告。
//
// 异常：
//  - 在无法创建工作目录、无法下载/拷贝输入、或读取 FASTA 失败时，会抛出 std::runtime_error。
//
// 输出（副作用）：
//  - 在 workdir/data/raw 中保存原始输入副本（或下载得到的文件）；
//  - 在 workdir/data/clean 中写入清洗后的 FASTA 文件和选中的共识候选（文件名见 config.hpp）；
//
// 约定：
//  - 该函数会尽量就地修改和写出数据以减少内存峰值；Top-K 选择器会保留 K 条完整记录在内存中。
uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n = 1000);


// ==============================================================
// alignConsensusSequence
//
// 说明：对外暴露的工具函数，用于对 `input_file` 中的未对齐共识序列执行多序列比对（MSA），
// 将比对结果写入 `output_file`。该函数通常在 `preprocessInputFasta` 之后调用，处理流程为：
//  1) 使用 `msa_cmd` 模板构造命令（模板可包含 {input} {output} {thread} 占位符）；
//  2) 在指定的工作目录 `workdir` 下执行该命令（通过 shell 或 cmd 模块），并等待其完成；
//  3) 函数会记录运行时间并对结果文件做基本检查（存在性与大小）。
//
// 参数：
//  - input_file: 要对齐的未对齐 FASTA（FilePath）
//  - output_file: 将写入比对结果的文件路径（FilePath）
//  - msa_cmd: 多序列比对命令模板字符串（例如 "mafft --auto {input} > {output}"）
//  - workdir: 用于运行命令时的当前工作目录（命令中的相对路径以此为基准）
//  - threads: 分配给 MSA 命令的线程数（传递给模板中的 {thread}），具体生效与否取决于所用 MSA 工具
//
// 性能提示：
//  - MSA 通常是 CPU 密集型、内存敏感的步骤，请根据目标机器调整 `threads` 和 MSA 工具的参数；
//  - 若 MSA 工具支持流式接口，可考虑在 future 中改为流式管道以减少磁盘 I/O。
//
// ==============================================================
void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, int threads);

// ==============================================================
// validateRefAlignedConsistency
//
// 说明：验证两个 FASTA 文件的序列一致性（删除 gap 后）。
// 
// 用途：当用户同时提供 -r/--ref（参考 FASTA）和 --ref-align（预对齐 MSA）时，
// 需要验证两个文件的序列内容是否匹配（允许序列顺序不同）。
//
// 设计：
// - 将 ref_fasta 中的所有序列按 ID 索引存储在 hash map 中（去 gap 后的版本）
// - 逐条读取 ref_aligned 中的序列，从 map 中查找匹配的序列并比较
// - 这样支持两个文件序列顺序不同的情况
//
// 参数：
//  - ref_fasta: 参考 FASTA 文件路径（-r/--ref）
//  - ref_aligned: 预对齐 MSA 文件路径（--ref-align）
//
// 异常：
//  - 序列 ID 不匹配：ref_aligned 中存在 ref_fasta 中没有的序列
//  - 序列内容不匹配：删除 gap 后的序列内容不相同
//  - 序列数不匹配：两个文件的序列总数不同
//  - 任何不匹配情况都会抛出 std::runtime_error
//
// ==============================================================
void validateRefAlignedConsistency(const FilePath& ref_fasta, const FilePath& ref_aligned);

#endif //HALIGN4_PREPROCESS_H
