#include "core/sampler.h"

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

    std::unique_ptr<Sampler> Clone(int seed);
};