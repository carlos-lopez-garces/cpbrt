#include "geometry.h"
#include "rng.h"

void StratifiedSample1D(Float *sample, int nSamples, RNG &rng, bool jitter);

void StratifiedSample2D(Point2f *sample, int xSamples, int ySamples, RNG &rng, bool jitter);

// Randomly permutes an array of samples.
template <typename T> void Shuffle(T *sample, int nSamples, int nDimensions, RNG &rng) {
    for (int i = 0; i < nSamples; ++i) {
        int other = i + rng.UniformUInt32(nSamples - i);

        for (int j = 0; j < nDimensions; ++j) {
            std::swap(
                sample[nDimensions * i + j],
                sample[nDimensions * other + j]
            );
        }
    }
}

void LatinHypercube(Float *sample, int nSamples, int nDimensions, RNG &rng);