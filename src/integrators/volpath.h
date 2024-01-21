#ifndef CPBRT_INTEGRATORS_VOLPATH_H
#define CPBRT_INTEGRATORS_VOLPATH_H

#include "core/cpbrt.h"
#include "core/camera.h"
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
        Film *film,
        const Bounds2i &pixelBounds,
        Float rrThreshold = 1,
        const std::string &lightSampleStrategy = "uniform"
    ) : SamplerIntegrator(camera, sampler, film),
        maxDepth(maxDepth),
        rrThreshold(rrThreshold),
        lightSampleStrategy(lightSampleStrategy)
    {}

    void Preprocess(const Scene &scene, Sampler &sampler);

    Spectrum Li(
        const RayDifferential &ray,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        int depth
    ) const;
};

VolPathIntegrator *CreateVolPathIntegrator(
    const ParamSet &params, 
    std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera,
    Film *film
);

#endif // CPBRT_INTEGRATORS_VOLPATH_H