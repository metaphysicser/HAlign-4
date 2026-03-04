#include "align.h"

#include <stdexcept>
#include <cctype>
#include <limits>
#include <string>

// CIGAR 编码/解码/转换与序列投影
// - CIGAR 压缩形式：unit=(len<<4)|op，编码 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
// - 关键操作：字符串转换、编码解码、插入/删除检测、序列投影对齐
// - 性能优化：避免 O(N^2)、预分配、反向填充

namespace cigar
{
    // 操作符与编码相互转换
    static inline uint32_t opCharToCode(const char op)
    {
        switch (op) {
        case 'M': return 0;
        case 'I': return 1;
        case 'D': return 2;
        case 'N': return 3;
        case 'S': return 4;
        case 'H': return 5;
        case 'P': return 6;
        case '=': return 7;
        case 'X': return 8;
        default:
            throw std::runtime_error(std::string("Unknown CIGAR op char: ") + op);
        }
    }

    static inline char opCodeToChar(const uint32_t code)
    {
        switch (code) {
        case 0: return 'M';
        case 1: return 'I';
        case 2: return 'D';
        case 3: return 'N';
        case 4: return 'S';
        case 5: return 'H';
        case 6: return 'P';
        case 7: return '=';
        case 8: return 'X';
        default:
            throw std::runtime_error("Unknown CIGAR op code: " + std::to_string(code));
        }
    }

    // 编码 (operation,len) -> CigarUnit：(len<<4)|op
    CigarUnit cigarToInt(char operation, uint32_t len)
    {
        constexpr uint32_t kMaxLen = (1u << 28) - 1u;
        if (len == 0 || len > kMaxLen) {
            throw std::runtime_error("cigarToInt: invalid length=" + std::to_string(len));
        }
        const uint32_t op = opCharToCode(operation);
        return (len << 4) | (op & 0x0Fu);
    }

    // 解码 CigarUnit -> (operation,len)
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
    {
        const uint32_t op = (cigar & 0x0Fu);
        len = (cigar >> 4);
        operation = opCodeToChar(op);
    }

    // 检测是否存在插入操作 'I'（O(n)，有短路优化）
    bool hasInsertion(const Cigar_t& cigar)
    {
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len;
            intToCigar(cu, op_char, len);
            if (op_char == 'I') return true;
        }
        return false;
    }

    // 压缩 Cigar_t -> SAM CIGAR string（如"10M5I3D"）
    std::string cigarToString(const Cigar_t& cigar)
    {
        if (cigar.empty()) return "";

        std::string out;
        out.reserve(cigar.size() * 5);  // 预分配避免扩容

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);
            out.append(std::to_string(len));
            out.push_back(op_char);
        }
        return out;
    }

    // SAM CIGAR string -> 压缩 Cigar_t（如"10M5I3D"解析为压缩向量）
    Cigar_t stringToCigar(const std::string& cigar_str)
    {
        Cigar_t result;
        if (cigar_str.empty() || cigar_str == "*") return result;

        result.reserve(cigar_str.size() / 2 + 1);

        uint64_t len_acc = 0;
        bool has_number = false;

        for (std::size_t i = 0; i < cigar_str.size(); ++i) {
            const unsigned char c = static_cast<unsigned char>(cigar_str[i]);

            if (std::isspace(c)) {
                continue;  // 跳过空白字符
            }

            if (std::isdigit(c)) {
                has_number = true;
                len_acc = len_acc * 10 + (c - '0');

                // 溢出保护
                if (len_acc > ((1ull << 28) - 1ull)) {
                    throw std::runtime_error("stringToCigar: op length overflow in '" + cigar_str + "'");
                }
                continue;
            }

            // 走到这里：c 不是数字也不是空白，认为是 op 字符
            if (!has_number || len_acc == 0) {
                throw std::runtime_error("stringToCigar: missing/invalid length before op in '" + cigar_str + "'");
            }

            const uint32_t len = static_cast<uint32_t>(len_acc);
            result.push_back(cigarToInt(static_cast<char>(c), len));

            // 重置状态准备下一个 op
            len_acc = 0;
            has_number = false;
        }

        // 字符串以数字结尾（缺少 op）属于格式错误
        if (has_number) {
            throw std::runtime_error("stringToCigar: trailing number without op in '" + cigar_str + "'");
        }

        return result;
    }

    // query 投影到 ref 坐标系：插入 gap '-'（从后往前填充避免 O(N^2)）
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar)
    {
        if (cigar.empty()) {
            return;
        }

        // 统计结果长度和消耗的 query 长度
        std::size_t out_len = 0;
        std::size_t consume_query = 0;

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (op_char == 'D' || op_char == 'N') {
                out_len += len;
            } else if (op_char == 'H' || op_char == 'P') {
                continue;
            } else {
                out_len += len;
                consume_query += len;
            }
        }

        assert(consume_query == query.size());

        // 备份原 query，分配输出空间
        std::string old = std::move(query);
        query.assign(out_len, '-');

        // 从后往前填充
        std::size_t w = out_len;
        std::size_t r = old.size();

        for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
            char op_char;
            uint32_t len = 0;
            intToCigar(*it, op_char, len);

            if (op_char == 'D' || op_char == 'N') {
                for (uint32_t i = 0; i < len; ++i) {
                    query[--w] = '-';
                }
                continue;
            }

            if (op_char == 'H' || op_char == 'P') {
                continue;
            }

            for (uint32_t i = 0; i < len; ++i) {
                assert(r > 0);
                query[--w] = old[--r];
            }
        }

        assert(w == 0);
        assert(r == 0);
    }

    // 按 CIGAR 调整 query：删除 I 操作碱基，为 D 操作添加 gap
    void delQueryToRefByCigar(std::string& query, const Cigar_t& cigar)
    {
        if (cigar.empty()) {
            return;
        }

        // 统计结果长度和 query 消耗长度
        std::size_t out_len = 0;
        std::size_t consume_query = 0;

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (op_char == 'I') {
                consume_query += len;
            } else if (op_char == 'H' || op_char == 'P') {
                continue;
            } else if (op_char == 'D' || op_char == 'N') {
                out_len += len;
            } else {
                out_len += len;
                consume_query += len;
            }
        }

        assert(consume_query == query.size());

        // 快速路径：无 I/D/N 操作则直接返回
        if (out_len == query.size()) {
            bool has_indel = false;
            for (const CigarUnit cu : cigar) {
                char op_char;
                uint32_t len = 0;
                intToCigar(cu, op_char, len);
                if (op_char == 'I' || op_char == 'D' || op_char == 'N') {
                    has_indel = true;
                    break;
                }
            }
            if (!has_indel) return;
        }

        // 备份原 query，分配输出空间
        std::string old = std::move(query);
        query.assign(out_len, '-');

        // 从后往前填充
        std::size_t w = out_len;
        std::size_t r = old.size();

        for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
            char op_char;
            uint32_t len = 0;
            intToCigar(*it, op_char, len);

            if (op_char == 'I') {
                // 跳过 I 操作的字符
                for (uint32_t i = 0; i < len; ++i) {
                    assert(r > 0);
                    --r;
                }
                continue;
            }

            if (op_char == 'H' || op_char == 'P') {
                continue;
            }

            if (op_char == 'D' || op_char == 'N') {
                // D/N 对应 gap（已初始化为 '-'）
                assert(w >= len);
                w -= len;
                continue;
            }

            // M/S/=/X：拷贝字符
            for (uint32_t i = 0; i < len; ++i) {
                assert(r > 0);
                assert(w > 0);
                query[--w] = old[--r];
            }
        }

        assert(w == 0);
        assert(r == 0);
    }

    // 智能追加 CIGAR 并合并相邻同类型操作
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add)
    {
        for (const CigarUnit cu : cigar_to_add) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (len == 0) continue;

            if (result.empty()) {
                result.push_back(cu);
            } else {
                char last_op_char;
                uint32_t last_len = 0;
                intToCigar(result.back(), last_op_char, last_len);

                if (last_op_char == op_char) {
                    constexpr uint32_t kMaxLen = (1u << 28) - 1u;
                    if (static_cast<uint64_t>(last_len) + len > kMaxLen) {
                        throw std::runtime_error("appendCigar: merged length overflow");
                    }
                    result.back() = cigarToInt(op_char, last_len + len);
                } else {
                    result.push_back(cu);
                }
            }
        }
    }

    // 计算 CIGAR 消耗的参考序列长度
    std::size_t getRefLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // M/D/N/=/X 消耗 ref
            if (op_char == 'M' || op_char == 'D' || op_char == 'N' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

    // 计算 CIGAR 消耗的查询序列长度
    std::size_t getQueryLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // M/I/S/=/X 消耗 query
            if (op_char == 'M' || op_char == 'I' || op_char == 'S' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

} // namespace cigar

