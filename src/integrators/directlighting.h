#ifndef CPBRT_INTEGRATORS_DIRECT_LIGHTING_H
#define CPBRT_INTEGRATORS_DIRECT_LIGHTING_H

#include <memory>

#include "core/cpbrt.h"
#include "core/integrator.h"
#include "core/scene.h"

// The strategy to choose depends on the number of samples taken per pixel: if many,
// then UniformSampleOne suffices (because many such single-sample, single-light
// samples will be taken anyway); if few, then UniformSampleAll to compensate.
enum class LightStrategy {
    // Take a number of samples from every light source and sum them.
    UniformSampleAll,
    // Take a single sample from only one light source chosen uniformly at random.
    UniformSampleOne
};

class DirectLightingIntegrator : public SamplerIntegrator {
private:
    const LightStrategy strategy;

    // Maximum recursion depth for rays traced for specular reflection or transmission.
    const int maxDepth;

    // Number of samples that are actually taken from the corresponding Light. Light::nSamples
    // tells how many should be taken from it, but the Sampler may round it up to a number
    // that better suits the sampling method. (Ordering follows Scene::lights's.)
    std::vector<int> nLightSamples;

public:
    DirectLightingIntegrator(
        LightStrategy strategy,
        int maxDepth,
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler,
        Film *film
    ) : SamplerIntegrator(camera, sampler, film),
        strategy(strategy),
        maxDepth(maxDepth) {}

    void Preprocess(const Scene &scene, Sampler &sampler);

    // In general, the estimate of incident radiance that the direct lighting integrator computes
    // uses samples from light sources that illuminate the intersection point directly. But if the
    // surface has a (perfect) specular BRDF, the light that bounces off of other surfaces coming
    // in the (perfect) specular reflection direction must be sampled as well: that's how the
    // intersected surface can reflect the image of nearby objects.
    Spectrum Li(
        const RayDifferential &ray,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        // The current ray tracing recursion depth for specular reflection and transmission.
        // If depth reaches maxDepth, the tracing recursion stops.
        int depth
    ) const;
};

DirectLightingIntegrator *CreateDirectLightingIntegrator(
    const ParamSet &params, 
    std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera,
    Film *film
);

#endif // CPBRT_INTEGRATORS_DIRECT_LIGHTING_H