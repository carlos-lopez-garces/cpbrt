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