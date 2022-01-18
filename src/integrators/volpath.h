#ifndef CPBRT_INTEGRATORS_VOLPATH_H
#define CPBRT_INTEGRATORS_VOLPATH_H

#include "core/cpbrt.h"
#include "core/integrator.h"
#include "core/lightdistribution.h"

class VolPathIntegrator : public SamplerIntegrator {
private:
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;

public:
    VolPathIntegrator(
        int maxDepth,
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler,
        const Bounds2i &pixelBounds,
        Float rrThreshold = 1,
        const std::string &lightSampleStrategy = "spatial"
    ) : SamplerIntegrator(camera, sampler, pixelBounds),
        maxDepth(maxDepth),
        rrThreshold(rrThreshold),
        lightSampleStrategy(lightSampleStrategy)
    {}

    void Preprocess(const Scene &scene, Sampler &sampler);
};

#endif // CPBRT_INTEGRATORS_VOLPATH_H