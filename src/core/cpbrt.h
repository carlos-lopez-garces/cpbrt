#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "error.h"

// Global forward declarations.

// Choose between float and double for most floating-point declarations
// at compile time.
#ifdef CPBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // CPBRT_FLOAT_AS_DOUBLE

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