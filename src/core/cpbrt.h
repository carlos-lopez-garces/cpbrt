#ifndef CPBRT_CORE_PBRT_H
#define CPBRT_CORE_PBRT_H

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

#if defined(_WIN32) || defined(_WIN64)
  #define CPBRT_IS_WINDOWS
#endif

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
    // Number of worker threads for parallel for loops.
    int nThreads = 0;
    bool quickRender = false;
    bool quiet = false;
    bool verbose = false;
    std::string imageFile;
    // Used by API state verification macros.
    bool cat = false;
    // Used by API state verification macros.
    bool toPly = false;
    // x0, x1, y0, y1.
    Float cropWindow[2][2];

    Options() {
        cropWindow[0][0] = 0;
        cropWindow[0][1] = 1;
        cropWindow[1][0] = 0;
        cropWindow[1][1] = 1;
    }
    
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

class Scene;
class Integrator;
class SamplerIntegrator;
template <typename T> class Vector2;
template <typename T> class Vector3;
template <typename T> class Point3;
template <typename T> class Point2;
template <typename T> class Normal3;
class Ray;
class RayDifferential;
template <typename T> class Bounds2;
template <typename T> class Bounds3;
class Transform;
class Interaction;
class SurfaceInteraction;
class Shape;
class Primitive;
class GeometricPrimitive;
class TransformedPrimitive;
template <int nSpectrumSamples> class CoefficientSpectrum;
class Camera;
struct CameraSample;
class ProjectiveCamera;
class Sampler;
class Filter;
class Film;
class FilmTile;
class BxDF;
class BRDF;
class BTDF;
class BSDF;
class Material;
template <typename T> class Texture;
class Medium;
class MediumInteraction;
struct MediumInterface;
class BSSRDF;
class SeparableBSSRDF;
class TabulatedBSSRDF;
struct BSSRDFTable;
class Light;
class VisibilityTester;
class AreaLight;
struct Distribution1D;
class Distribution2D;
class RNG;
class ProgressReporter;
class MemoryArena;
template <typename T, int logBlockSize = 2> class BlockedArray;
class Matrix4x4;
class ParamSet;
template <typename T>
struct ParamSetItem;
class TextureParams;

// Global constants.

// Only so that you don't have to type std::numeric_limits<Float> ... 
#ifdef _MSC_VER
#define Infinity std::numeric_limits<Float>::infinity()
#define MaxFloat std::numeric_limits<Float>::max()
#else
static CPBRT_CONSTEXPR Float Infinity = std::numeric_limits<Float>::infinity();
static CPBRT_CONSTEXPR Float MaxFloat = std::numeric_limits<Float>::max();
#endif // _MSC_VER

// The machine epsilon is half the magnitude of the minimum spacing between real numbers
// that are representable in floating-point. The spacing between consecutive floating-
// point numbers varies with the range of powers of 2: the spacing between numbers in 
// the [2^ek, 2^ek+1] range is smaller than the spacing in the [2^ek+1, 2^ek+2] range.
// The range closest to 0 is where the spacing is the smallest.
//
// Half the epsilon is the upper bound on the floating-point rounding error of an
// arithmetic operation, when the result falls in the power of 2 range closest to 0. An
// arithmetic operation rounds the resulting real number up or down to the closest
// floating-point number; the difference between the real number and the closest floating-
// point number is the rounding error and it can't be larger than half the epsilon (again,
// when the real number falls in the power of 2 range closest to 0; for other ranges,
// the upper bound of the rounding error scales in proportion to the spacing of that range). 
#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)
#else
static PBRT_CONSTEXPR Float MachineEpsilon =
    std::numeric_limits<Float>::epsilon() * 0.5;
#endif

static const Float Pi      = 3.14159265358979323846;
static const Float InvPi   = 0.31830988618379067154;
static const Float Inv2Pi  = 0.15915494309189533577;
static const Float Inv4Pi  = 0.07957747154594766788;
static const Float PiOver2 = 1.57079632679489661923;
static const Float PiOver4 = 0.78539816339744830961;
static const Float Sqrt2   = 1.41421356237309504880;
static CPBRT_CONSTEXPR Float ShadowEpsilon = 0.0001f;

// Dynamic memory allocation on the stack.
#ifdef _MSC_VER
#define alloca _alloca
#endif

#ifndef CPBRT_L1_CACHE_LINE_SIZE
#define CPBRT_L1_CACHE_LINE_SIZE 64
#endif

// Global inline functions.

inline Float GammaCorrect(Float value) {
    if (value <= 0.0031308f) return 12.92f * value;
    return 1.055f * std::pow(value, (Float)(1.f / 2.4f)) - 0.055f;
}

inline Float InverseGammaCorrect(Float value) {
    if (value <= 0.04045f) return value * 1.f / 12.92f;
    return std::pow((value + 0.055f) * 1.f / 1.055f, (Float)2.4f);
}

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

inline float NextFloatUp(float v) {
    // There's no float higher than infinity.
    if (std::isinf(v) && v > 0.) return v;

    // -0 and 0 are different in their floating-point representation: -0 has a 1 in the sign bit.
    // Since they differ only in sign, their floating-point representations are not consecutive,
    // so incrementing -0 wouldn't yield 0.
    if (v == -0.f) v = 0.f;

    // Convert to a base 2 integer so that the increment and decrement operators yield the next
    // float.
    uint32_t ui = FloatToBits(v);
    if (v >= 0) {
        ++ui;
    } else {
        --ui;
    }
    return BitsToFloat(ui);
}

inline float NextFloatDown(float v) {
    // There's no float lower than -infinity.
    if (std::isinf(v) && v < 0.) return v;

    // -0 and 0 are different in their floating-point representation: -0 has a 1 in the sign bit.
    // Since they differ only in sign, their floating-point representations are not consecutive,
    // so decrementing 0 wouldn't yield -0.
    if (v == 0.f) v = -0.f;

    // Convert to a base 2 integer so that the increment and decrement operators yield the next
    // float.
    uint32_t ui = FloatToBits(v);
    if (v > 0) {
        --ui;
    } else {
        ++ui;
    }

    return BitsToFloat(ui);
}

inline double NextFloatUp(double v, int delta = 1) {
    if (std::isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;
    uint64_t ui = FloatToBits(v);
    if (v >= 0.)
        ui += delta;
    else
        ui -= delta;
    return BitsToFloat(ui);
}

inline double NextFloatDown(double v, int delta = 1) {
    if (std::isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint64_t ui = FloatToBits(v);
    if (v > 0.)
        ui -= delta;
    else
        ui += delta;
    return BitsToFloat(ui);
}

inline Float gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

// TODO.
inline Float Erf(Float x) {
    // constants
    Float a1 = 0.254829592f;
    Float a2 = -0.284496736f;
    Float a3 = 1.421413741f;
    Float a4 = -1.453152027f;
    Float a5 = 1.061405429f;
    Float p = 0.3275911f;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    Float t = 1 / (1 + p * x);
    Float y =
        1 -
        (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

// TODO.
inline Float ErfInv(Float x) {
    Float w, p;
    x = Clamp(x, -.99999f, .99999f);
    w = -std::log((1 - x) * (1 + x));
    if (w < 5) {
        w = w - 2.5f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    } else {
        w = std::sqrt(w) - 3;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

template <typename T> inline bool IsPowerOf2(T v) {
    // A power of 2 in binary notation is a single 1 followed by 0s.
    // Let v = 4 = 100b. Then v-1 = 011b and v & (v-1) = 100b & 011b = 0.
    return v && !(v & (v-1));
}

extern Options CpbrtOptions;

#endif // CPBRT_CORE_PBRT_H