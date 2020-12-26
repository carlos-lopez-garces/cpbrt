#include <algorithm>
#ifdef CPBRT_HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#ifdef CPBRT_HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <memory>
#include <string>
#include <vector>

#include "error.h"

// Global macros.

// Dynamic memory allocation on the stack.
#define ALLOCA(TYPE, COUNT) (TYPE *) alloca((COUNT) * sizeof(TYPE))

// Global forward declarations.

// Choose between float and double for most floating-point declarations
// at compile time.
#ifdef CPBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // CPBRT_FLOAT_AS_DOUBLE

struct Options {

};

// Rather than being an abstract data type, Spectrum is an alias of the chosen
// concrete implementation. This allows 1) other types and structs to have members
// of type Spectrum that aren't pointers to dynamically-allocated objects, and
// 2) for inlining short functions, which would otherwise be virtual invocations
// resolved at run-time.
class RGBSpectrum;
class SampledSpectrum;
#ifdef CPBRT_SAMPLED_SPECTRUM
typedef SampledSpectrum Spectrum;
#else
typedef RGBSpectrum Spectrum;
#endif // CPBRT_SAMPLED_SPECTRUM

// Global constants.

// Only so that you don't have to type std::numeric_limits<Float> ... 
#ifdef _MSC_VER
#define Infinity std::numeric_limits<Float>::infinity()
#else
static CPBRT_CONSTEXPR Float Infinity = std::numeric_limits<Float>::infinity();
#endif // _MSC_VER

static /*CPBRT_CONSTEXPR*/ const Float Pi = 3.14159265358979323846;

// Dynamic memory allocation on the stack.
#ifdef _MSC_VER
#define alloca _alloca
#endif

// Global inline functions.

inline Float Radians(Float deg) {
    return deg * Pi / 180.f;
}

inline Float Degrees(Float rad) {
    return rad * 180.f / Pi;
}

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
    if (val < low) {
        return low;
    }
    else if (val > high) {
        return high;
    }
    return val;
}

inline Float Lerp(Float t, Float v1, Float v2) {
    return (1-t)*v1 + t*v2;
}

// Given the size of an interval and a predicate, FindInterval partitions the [0,size]
// interval of indices into n=size subintervals and computes the index of the left endpoint
// of the subinterval for which the input predicate is true.
//
// This partitioned interval of monotonically-increasing integer indices is imposed over
// the possibly irregularly-partitioned interval defined implicitly by the predicate. The
// subinterval size, as well as the value that corresponds to an endpoint index, are known
// only by the predicate.
template <typename Predicate>
int FindInterval(int size, const Predicate &pred) {
    int first = 0;
    int len = size;

    // Do a binary search guided by the predicate.
    while (len > 0) {
        int half = len >> 1;
        int middle = first + half;

        // Bisect interval based on value of predicate at middle.
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }

    int leftEndpoint = Clamp(first-1, 0, size-2);
    return leftEndpoint;
}

inline uint32_t FloatToBits(float f) {
    uint32_t bits;
    memcpy(&bits, &f, sizeof(float));
    return bits;
}

inline uint64_t FloatToBits(double f) {
    uint64_t bits;
    memcpy(&bits, &f, sizeof(double));
    return bits;
}

inline float BitsToFloat(uint32_t bits) {
    float f;
    memcpy(&f, &bits, sizeof(uint32_t));
    return f;
}

inline double BitsToFloat(uint64_t bits) {
    double f;
    memcpy(&f, &bits, sizeof(uint64_t));
    return f;
}