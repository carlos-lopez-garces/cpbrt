#include <atomic>

#include "cpbrt.h"

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