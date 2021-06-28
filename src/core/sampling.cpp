#include "sampling.h"
#include "geometry.h"
#include "shape.h"

void StratifiedSample1D(Float *sample, int nSamples, RNG &rng, bool jitter) {
    // Normalized size of stratum.
    Float stratumScale = (Float) 1 / nSamples;

    for (int i = 0; i < nSamples; ++i) {
        // Offset within stratum: 0.5 for regular, uniform stratified sampling; a random
        // offset from [0,1) for jittered stratified sampling.
        Float offset = jitter ? rng.UniformFloat() : 0.5f;

        // Draw the ith sample from the ith stratum.
        // 
        // Each stratum is a subinterval of [0,1). The start point of the ith stratum is
        // i*stratumScale; the midpoint is (i + 0.5)*stratumScale. If jittered, the offset
        // may choose any point within the stratum's subinterval.
        sample[i] = std::min((i + offset) * stratumScale, OneMinusEpsilon);
    }
}

void StratifiedSample2D(Point2f *sample, int xSamples, int ySamples, RNG &rng, bool jitter) {
    Float xStratumScale = (Float) 1 / xSamples;
    Float yStratumScale = (Float) 1 / ySamples;

    for (int y = 0; y < ySamples; ++y) {
        for (int x = 0; x < xSamples; ++x) {
            // Offsets within stratum along each direction are independent and random.
            Float xOffset = jitter ? rng.UniformFloat() : 0.5f;
            Float yOffset = jitter ? rng.UniformFloat() : 0.5f;
            sample->x = std::min((x + xOffset) * xStratumScale, OneMinusEpsilon);
            sample->y = std::min((y + yOffset) * yStratumScale, OneMinusEpsilon);
            ++sample;
        }
    }
}

void LatinHypercube(Float *sample, int nSamples, int nDimensions, RNG &rng) {
    // Generate samples along diagonal.
    Float reciprocalNSamples = (Float) 1 / nSamples;
    for (int i = 0; i < nSamples; ++i) {
        for (int j = 0; j < nDimensions; ++j) {
            Float jitteredSample = (i + rng.UniformFloat()) * reciprocalNSamples;
            sample[nDimensions * i + j] = std::min(jitteredSample, OneMinusEpsilon);
        }
    }

    // Permute samples in each dimension.
    for (int i = 0; i < nDimensions; ++i) {
        for (int j = 0; j < nSamples; ++j) {
            int other = j + rng.UniformUInt32(nSamples - j);
            std::swap(sample[nDimensions * j + i], sample[nDimensions * other  + i]);
        }
    }
}

// Samples the unit disk uniformly using a pair of [0,1] uniform random numbers.
// The output sample (X,Y) is the cartesian coordinate of a point on the disk, obtained
// by transforming the polar coordinate of a sample of a 2D random variable (r, theta).
//
// The value of r comes from sampling the marginal PDF of r, which is accomplished by the
// inversion method by evaluating the inverse of the CDF of r with the value of a uniform
// random variable. (The CDF of r is obtained by integrating the marginal PDF.) 
//
// The value of theta comes from sampling the conditional density function p(theta|r) given
// the value of r obtained before. p(theta|r) is sampled by the inversion method by evaluating
// the inverse with the value of a uniform random variable.
Point2f UniformSampleDisk(const Point2f &u) {
    Float r = std::sqrt(u[0]);
    Float theta = 2 * Pi * u[1];
    // Map the polar coordinate to a cartesian coordinate.
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

// Samples the unit disk using a concentric mapping of the unit square, which transforms a
// uniformly distributed random point on the unit square to a point on the unit disk. 
Point2f ConcentricSampleDisk(const Point2f &u) {
    // Map uniform random numbers to the unit square [-1,1]^2.
    Point2f uOffset = 2.f * u - Vector2f(1, 1);

    // Handle degeneracy at the origin.
    if (uOffset.x == 0 && uOffset.y == 0) {
        return Point2f(0, 0);
    }

    // Apply concentric mapping from the unit square to the unit disk. This mapping turns
    // wedges of the square into slices of the disk.
    Float theta, r;
    if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
        r = uOffset.x;
        theta = PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
    }

    // Map the polar coordinate to a cartesian coordinate.
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

Vector3f UniformSampleHemisphere(const Point2f &u) {
    Float z = u[0];
    Float r = std::sqrt(std::max((Float)0, (Float)1. - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformHemispherePdf() {
    return Inv2Pi;
}

// TODO: explain.
Vector3f UniformSampleSphere(const Point2f &u) {
    Float z = 1 - 2 * u[0];
    Float r = std::sqrt(std::max((Float)0, (Float)1 - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

// TODO: explain.
Float UniformSpherePdf() {
    return Inv4Pi;
}

// TODO: explain.
Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax) {
    Float cosTheta = ((Float)1 - u[0]) + u[0] * cosThetaMax;
    Float sinTheta = std::sqrt((Float)1 - cosTheta * cosTheta);
    Float phi = u[1] * 2 * Pi;
    return Vector3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
}

// TODO: explain.
Vector3f UniformSampleCone(
    const Point2f &u,
    Float cosThetaMax,
    const Vector3f &x,
    const Vector3f &y,
    const Vector3f &z
) {
    Float cosTheta = Lerp(u[0], cosThetaMax, 1.f);
    Float sinTheta = std::sqrt((Float)1. - cosTheta * cosTheta);
    Float phi = u[1] * 2 * Pi;
    return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y + cosTheta * z;
}

// TODO: explain.
Float UniformConePdf(Float cosThetaMax) {
    return 1 / (2 * Pi * (1 - cosThetaMax));
}