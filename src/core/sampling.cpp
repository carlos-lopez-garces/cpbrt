#include "sampling.h"

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