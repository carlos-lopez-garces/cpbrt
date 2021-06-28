#ifndef CPBRT_INTEGRATORS_PATH_H
#define CPBRT_INTEGRATORS_PATH_H

#include "core/cpbrt.h"
#include "core/integrator.h"

class PathIntegrator : public SamplerIntegrator {
private:
    // Constant-length path sampling and path estimate evaluation termination condition.
    // Russian roulette might terminate the path earlier, though.
    const int maxDepth;

public:
    PathIntegrator(
        int maxDepth,
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler
    ) : SamplerIntegrator(camera, sampler),
        maxDepth(maxDepth)
    {}

    Spectrum Li(
        const RayDifferential &r,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        int depth
    ) const;
};

PathIntegrator *CreatePathIntegrator(
    const ParamSet &params,
    std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera
);

#endif // CPBRT_INTEGRATORS_PATH_H