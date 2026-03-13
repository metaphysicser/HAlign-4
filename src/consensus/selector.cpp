#include "consensus.h"

// ---------------- TopKLongestSelector 实现说明（详细中文注释） ----------------
// TopKLongestSelector 的职责：在一次线性扫描（streaming）中，选出长度最长的前 K 条序列，
// 并保持在长度相同情况下的稳定性（优先保留更早出现的序列）。
//
// 设计目标与动机：
// - 流式（streaming）处理：不需要将所有序列都加载到内存；适用于 N 很大（百万级或更高）的输入。
// - 空间受限：只使用 O(K) 的额外空间来保存候选；K 通常远小于 N（例如 K 为几百或几千）。
// - 时间高效：每读到一条序列，做 O(log K) 的堆操作；总体 O(N log K) 时间复杂度。
// - 稳定性：当序列长度相同时，保留先出现的那条（稳定的选择有助于可复现性和简单调试）。
//
// 数据结构与比较准则：
// - 使用最小堆（min-heap）：堆顶保存当前 K 个候选中的“最差”元素（即最短或在同样长度下最晚出现），
//   这样新的更好的候选只需与堆顶比较，若更好则替换堆顶并下沉，保持堆的大小为 K。
// - Item 包含字段：len（序列长度）、order（输入顺序计数器，用于稳定性比较）、rec（SeqRecord 本身）。
// - 比较准则：词典序 (len, -order) 或在实现中用两个函数表达：
//     * worseThan(a,b)：若 a 在排序上比 b 更“差”（即应该排在前面的更差元素放在堆顶），返回 true；
//       这里定义为：长度更短 -> 更差；若长度相等，order 更大（出现更晚） -> 更差。
//     * betterThan(a,b)：相反的判断，用于判定候选是否优于当前堆顶。
//
// 稳定性说明：
// - 我们通过 order_counter_ 递增记录读到每条记录的先后顺序。
// - 当长度相同时，出现顺序更早的记录被认为更好（order 值更小）。因此在相同长度时，先到的记录优先保留，
//   从而保证最终结果的稳定性（同样的输入顺序给出一致结果）。
//
// 内存与性能权衡：
// - 内存：保留 K 条 SeqRecord 的内存为 O(K * avg_len)。若 avg_len 很大（例如每条序列很长），并且 K 很大，
//   可能造成较高内存占用；若内存敏感，可以仅保留索引/元数据或做外部存储。
// - 性能：堆操作为 O(log K)。若 K 很小（例如 10、100），堆操作会非常快。若 K 挺大（接近 N），应考虑其它策略（如 partial sort 或外部排序）。
//
// 并发/线程安全：
// - 本实现不是线程安全的；TopKLongestSelector 假定在单线程上下文中被访问（例如在读文件的主线程中逐条调用 consider）。
// - 若希望多线程并行处理输入（多线程读序列并同时 consider），需要对 selector 加锁或采用线程本地局部 TopK 后再归并（推荐）：
//     * 每个线程维护自己的 TopK（thread-local），在处理完一段数据后，把每个线程的 TopK 合并到全局 TopK（合并成本 O(T * K log K)）。
//
// 可替代实现/扩展：
// - 若只需要最终 TopK 的 Unstable 版本（不关心稳定性），可以使用 std::priority_queue + 简化比较函数。
// - 若 K 很大且内存不足，可采用外部归并（external merge）或对输入进行两轮扫描（第一轮抽样估计阈值，第二轮过滤）。
// - 若需要支持按其他度量（例如按质量评分或某种复合度量），只需调整 Item 的比较函数即可。
//
// 边界与异常情况处理：
// - 若 k_ == 0：selector 在 consider 时直接返回，不保存任何记录；takeSortedDesc 返回空向量。
// - 如果输入中存在超长序列导致 len 值溢出到 size_t 的极端情况（非常不现实），代码会按 size_t 语义处理。
// - 本实现假定 seq_io::SeqRecord 的移动语义正确（即可以通过 std::move(rec) 将内存转移到堆中而不是复制）。


TopKLongestSelector::TopKLongestSelector(std::size_t k)
        : k_(k)
    {
        heap_.reserve(k_);
    }

    // 重置选择器为新的 K 值，清空内部状态
    void TopKLongestSelector::reset(std::size_t k)
    {
        k_ = k;
        order_counter_ = 0; // order_counter_ 用于记录元素出现的先后顺序，以便在长度相同时保持稳定性
        heap_.clear();
        heap_.reserve(k_);
    }

    // 当前堆中已保存的元素数量
    std::size_t TopKLongestSelector::size() const
    {
        return heap_.size();
    }

    // selector 的容量 K
    std::size_t TopKLongestSelector::capacity() const
    {
        return k_;
    }

    // 是否为空
    bool TopKLongestSelector::empty() const
    {
        return heap_.empty();
    }

    // 比较函数：判断 a 是否比 b 更 "差"
    // 语义：在最小堆中，“更差”的元素会被置于堆顶（即应当被替换或弹出）。
    bool TopKLongestSelector::worseThan(const Item& a, const Item& b)
    {
        if (a.len != b.len) return a.len < b.len;      // 长度更短的被认为更差
        return a.order > b.order;                      // 同长度时，后来出现的（order 值更大）被认为更差
    }

    // 判断候选项 cand 是否比当前最坏项 worst 更好（用于决定是否替换堆顶）
    bool TopKLongestSelector::betterThan(const Item& cand, const Item& worst)
    {
        if (cand.len != worst.len) return cand.len > worst.len; // 更长则更好
        // 同长度时，我们偏向保留更早出现的（即 order 更小的项优先）
        return cand.order < worst.order;
    }

    // 上浮操作：当新元素插入到堆尾时调用，若其比父节点更差则与父节点交换直到堆序恢复
    // 说明：这里实现的是一个 min-heap 的上浮（维护堆顶为最差元素），使得堆顶始终为当前最差
    void TopKLongestSelector::siftUp(std::size_t idx)
    {
        while (idx > 0) {
            const std::size_t parent = (idx - 1) / 2;
            // min-heap：如果当前更“差”（更小），就上浮（与父节点交换）
            if (worseThan(heap_[idx], heap_[parent])) {
                std::swap(heap_[idx], heap_[parent]);
                idx = parent;
            } else {
                break;
            }
        }
    }

    // 下沉操作：当堆顶被替换后调用，把新的堆顶下沉到合适位置
    // 通过比较左右孩子找到更"差"的孩子并与当前交换，直到堆序恢复
    void TopKLongestSelector::siftDown(std::size_t idx)
    {
        const std::size_t n = heap_.size();
        while (true) {
            const std::size_t left = idx * 2 + 1;
            if (left >= n) break; // 无孩子，停止

            const std::size_t right = left + 1;

            // 选择更“差”的孩子（min-heap 中更差的孩子应当成为交换对象）
            std::size_t worst_child = left;
            if (right < n && worseThan(heap_[right], heap_[left])) {
                worst_child = right;
            }

            // 如果孩子比当前更“差”，则下沉（交换）
            if (worseThan(heap_[worst_child], heap_[idx])) {
                std::swap(heap_[idx], heap_[worst_child]);
                idx = worst_child;
            } else {
                break;
            }
        }
    }

    // 考虑一条新记录：将其包装为 Item 并尝试加入堆中
    // 逻辑：
    // - 如果 k_ == 0，直接返回（不保存任何记录）
    // - 如果堆未满，则直接 push 并上浮
    // - 否则比较候选与堆顶（最差）；若候选更好，则替换堆顶并下沉
    //
    // 性能优化（关键）：
    // - 先检查长度和 order，再决定是否移动 SeqRecord
    // - 避免对 99.99% 的不合格记录进行无效的字符串移动构造
    // - 对于 100 万条记录，K=100 的场景，可减少 ~99.99 万次字符串移动
    void TopKLongestSelector::consider(seq_io::SeqRecord rec)
    {
        if (k_ == 0) return;

        const std::size_t cand_len = rec.seq.size();
        const std::uint64_t cand_order = order_counter_++;

        if ((double)rec.n_num / (double)rec.seq.size() > 0.01) {
            return;
        }

        // 堆未满：直接接受
        if (heap_.size() < k_) {
            Item item;
            item.len = cand_len;
            item.order = cand_order;
            item.rec = std::move(rec);
            heap_.push_back(std::move(item));
            siftUp(heap_.size() - 1);
            return;
        }

        // 堆已满：先用长度和 order 判断是否比堆顶更好
        // 关键优化：避免无效的 SeqRecord 移动
        const Item& worst = heap_[0];

        // 快速判断：候选是否比堆顶更好
        bool is_better;
        if (cand_len != worst.len) {
            is_better = (cand_len > worst.len);  // 更长则更好
        } else {
            is_better = (cand_order < worst.order);  // 同长度时，更早出现更好
        }

        // 只有确认更好时，才移动 SeqRecord 并替换堆顶
        if (is_better) {
            heap_[0].len = cand_len;
            heap_[0].order = cand_order;
            heap_[0].rec = std::move(rec);
            siftDown(0);
        }
        // 否则：rec 被销毁，但避免了不必要的移动到 Item 再销毁的开销
    }

    // 将保留的堆内容按长度降序输出（同长度按出现顺序升序，保证稳定性）
    // 返回值：vector<SeqRecord>，其中包含已选择的记录（按长度从长到短排序）
    std::vector<seq_io::SeqRecord> TopKLongestSelector::takeSortedDesc()
    {
        // 为避免复制大量序列数据，move 出 heap_ 的内容到一个临时 vector 上进行排序
        std::vector<Item> items = std::move(heap_);
        heap_.clear();
        heap_.shrink_to_fit();
        heap_.reserve(k_);

        // 自定义排序：长度降序；若长度相同则按原始出现顺序升序（更早的排在前面）
        std::sort(items.begin(), items.end(),
                  [](const Item& a, const Item& b) {
                      if (a.len != b.len) return a.len > b.len;
                      return a.order < b.order;
                  });

        // 将 SeqRecord 从 items 移出到返回向量中
        std::vector<seq_io::SeqRecord> out;
        out.reserve(items.size());
        for (auto& it : items) {
            out.push_back(std::move(it.rec));
        }
        return out;
    }
