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

// 预处理输入 FASTA，并返回处理的序列数量（total records processed）
//
// 参数：
//  - input_path: 输入 FASTA 文件路径，支持本地路径或远程 URL（字符串）。
//  - workdir:   工作目录（字符串），会在该目录下创建 data/raw 和 data/clean 子目录并写入中间文件。
//  - consensus_num:    要为后续共识选择的序列数量（Top-K，按序列长度挑选），默认为 1000。
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
uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int consensus_num = 1000);


// ==============================================================
// alignConsensusSequence
//
// 说明：对外暴露的工具函数，用于对 `input_file` 中的未对齐共识序列执行多序列比对（MSA），
// 将比对结果写入 `output_file`。该函数通常在 `preprocessInputFasta` 之后调用，处理流程为：
//  1) 使用 `msa_tool` 模板构造命令（模板可包含 {input} {output} {thread} 占位符）；
//  2) 在指定的工作目录 `workdir` 下执行该命令（通过 shell 或 cmd 模块），并等待其完成；
//  3) 函数会记录运行时间并对结果文件做基本检查（存在性与大小）。
//
// 参数：
//  - input_file: 要对齐的未对齐 FASTA（FilePath）
//  - output_file: 将写入比对结果的文件路径（FilePath）
//  - msa_tool: 多序列比对命令模板字符串（例如 "mafft --auto {input} > {output}"）
//  - workdir: 用于运行命令时的当前工作目录（命令中的相对路径以此为基准）
//  - threads: 分配给 MSA 命令的线程数（传递给模板中的 {thread}），具体生效与否取决于所用 MSA 工具
//  - verbose: 是否输出详细命令和文件检查日志；批量调用时可设为 false 以避免刷屏
//
// 性能提示：
//  - MSA 通常是 CPU 密集型、内存敏感的步骤，请根据目标机器调整 `threads` 和 MSA 工具的参数；
//  - 若 MSA 工具支持流式接口，可考虑在 future 中改为流式管道以减少磁盘 I/O。
//
// ==============================================================
void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_tool, int threads, bool verbose = true);

// ==============================================================
// validateRefAlignedConsistency
//
// 说明：验证两个 FASTA 文件的序列一致性（删除 gap 后）。
// 
// 用途：当用户同时提供 -r/--reference（参考 FASTA）和 --reference-msa（预对齐 MSA）时，
// 需要验证两个文件的序列内容是否匹配（允许序列顺序不同）。
//
// 设计：
// - 将 ref_fasta 中的所有序列按 ID 索引存储在 hash map 中（去 gap 后的版本）
// - 逐条读取 reference_msa 中的序列，从 map 中查找匹配的序列并比较
// - 这样支持两个文件序列顺序不同的情况
//
// 参数：
//  - ref_fasta: 参考 FASTA 文件路径（-r/--reference）
//  - reference_msa: 预对齐 MSA 文件路径（--reference-msa）
//
// 异常：
//  - 序列 ID 不匹配：reference_msa 中存在 ref_fasta 中没有的序列
//  - 序列内容不匹配：删除 gap 后的序列内容不相同
//  - 序列数不匹配：两个文件的序列总数不同
//  - 任何不匹配情况都会抛出 std::runtime_error
//
// ==============================================================
void validateRefAlignedConsistency(const FilePath& ref_fasta, const FilePath& reference_msa);

std::array<int8_t, 25> readScoreMatrixFile(const std::string& path);

#endif //HALIGN4_PREPROCESS_H
