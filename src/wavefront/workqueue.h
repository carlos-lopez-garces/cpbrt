#ifndef CPBRT_WAVEFRONT_WORKQUEUE_H
#define CPBRT_WAVEFRONT_WORKQUEUE_H

#include <atomic>

// WorkQueue generates an SOA (structure-of-arrays) data access pattern.
template <typename WorkItem>
class WorkQueue : public SOA<WorkItem> {
private:
    // Number of objects stored in the queue.
    std::atomic<int> size{0};
protected:
public:
    WorkQueue(int n, Allocator alloc) : SOA<WorkItem>(n, alloc) {}

    int Size() const {
        // memory_order_relaxed: there are no synchronization or ordering
        // constraints imposed on other reads or writes, only this operation's
        // atomicity is guaranteed.
        return size.load(std::memory_order_relaxed);
    }

    void Reset() {
        // memory_order_relaxed: there are no synchronization or ordering
        // constraints imposed on other reads or writes, only this operation's
        // atomicity is guaranteed.
        size.store(0, std::memory_order_relaxed);
    }
};

#endif // CPBRT_WAVEFRONT_WORKQUEUE_H