#ifndef CPBRT_CORE_PARALLEL_H
#define CPBRT_CORE_PARALLEL_H

#include <atomic>
#include <functional>

#include "cpbrt.h"
#include "geometry.h"

// TODO: std::atomic supports floating-point types since C++20.
class AtomicFloat {
private:
#ifdef CPBRT_FLOAT_AS_DOUBLE
    std::atomic<uint64_t> bits;
#else
    std::atomic<uint32_t> bits;
#endif

public:
    explicit AtomicFloat(Float v = 0) {
        bits = FloatToBits(v);
    }

    // Implicit conversion to Float.
    operator Float() const {
        return BitsToFloat(bits);
    }

    Float operator=(Float v) {
        bits = FloatToBits(v);
        return v;
    }

    void Add(Float v) {
#ifdef CPBRT_FLOAT_AS_DOUBLE
        uint64_t oldBits = bits;
        uint64_t newBits;
#else
        uint32_t oldBits = bits;
        uint32_t newBits;
#endif

        do {
            newBits = FloatToBits(BitsToFloat(oldBits) + v);

           // If bits == oldBits, set bits = newBits. Otherwise retry.
        } while (!bits.compare_exchange_weak(oldBits, newBits));
    }
};

// Thread-local variable that tells the thread what its identity is.
// External linkage so that other translation units can reference it.
extern CPBRT_THREAD_LOCAL int ThreadIndex;

void ParallelFor(
    // Wraps a callable: a lambda expression, a function, a pointer to a function, etc.
    // The callable receives the iteration's index.
    const std::function<void(int64_t)> &func,
    // Total number of iterations.
    int64_t  count,
    int chunkSize = 1
);

void ParallelFor2D(
    std::function<void(Point2i)> func,
    const Point2i &count
);

// Signals all threads to exit and cleans up resources.
void ParallelCleanup();

int NumSystemCores();

#endif // CPBRT_CORE_PARALLEL_H