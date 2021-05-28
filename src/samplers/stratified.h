#ifndef CPBRT_SAMPLERS_STRATIFIED_H
#define CPBRT_SAMPLERS_STRATIFIED_H

#include "core/sampler.h"
#include "core/rng.h"

class StratifiedSampler : PixelSampler {
private:
    const int xPixelSamples;
    const int yPixelSamples;
    const int jitterSamples;

public:
    StratifiedSampler(
        int xPixelSamples,
        int yPixelSamples,
        bool jitterSamples,
        int nSampledDimensions
    ) : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
        xPixelSamples(xPixelSamples),
        yPixelSamples(yPixelSamples),
        jitterSamples(jitterSamples) {}

    void StartPixel(const Point2i &p);

    std::unique_ptr<Sampler> Clone(int seed);
};

#endif // CPBRT_SAMPLERS_STRATIFIED_H