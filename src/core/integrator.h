#ifndef CPBRT_CORE_INTEGRATOR_H
#define CPBRT_CORE_INTEGRATOR_H

#include <memory>

#include "cpbrt.h"

// Evaluates the direct lighting outgoing radiance / scattering equation at the 
// intersection point with direct lighting contributions from all light sources in
// the scene, using Monte Carlo integration.
Spectrum UniformSampleAllLights(
    const Interaction &it,
    const Scene &scene,
    MemoryArena &arena,
    Sampler &sampler,
    // The number of samples to take by light source.
    const std::vector<int> &nLightSamples,
    // Whether to account for the effects of volumetric attenuation.
    bool handleMedia = false
);

// Evaluates the direct lighting outgoing radiance / scattering equation at the
// intersection point by taking a single sample from a single light source chosen
// uniformly at random. 
Spectrum UniformSampleOneLight(
    const Interaction &it,
    const Scene &scene,
    MemoryArena &arena,
    Sampler &sampler,
    // Whether to account for the effects of volumetric attenuation.
    bool handleMedia = false
);

// Computes a single Monte Carlo estimate of radiance emitted by the given light
// source and received by the surface or medium represented by the given Interaction,
// using multiple importance sampling.
Spectrum EstimateDirect(
    const Interaction &it,
    // Uniform random sample for sampling the BSDF.
    const Point2f &uScattering,
    const Light &light,
    // Uniform random sample for sampling the light source.
    const Point2f &uLight,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    // Whether to account for the effects of volumetric attenuation.
    bool handleMedia = false,
    // Whether to consider perfectly specular lobes.
    bool specular = false
);

class Integrator {
public:
    virtual void Render(const Scene &scene) = 0;
};

class SamplerIntegrator : public Integrator {
private:
    // Supplies sample points on the image at which to compute radiance.
    std::shared_ptr<Sampler> sampler;

protected:
    std::shared_ptr<const Camera> camera;

public:
    SamplerIntegrator(
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler
    ) : camera(camera), sampler(sampler) {}

    void Render(const Scene &scene);

    // Executes before the main rendering loop.
    virtual void Preprocess(const Scene &scene, Sampler &sampler) {}

    // Samples the incident radiance function using the given ray (and differentials).
    virtual Spectrum Li(
        const RayDifferential &ray,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        // Number of ray bounces from the camera that have occurred up until this call.
        int depth = 0
    ) const = 0;
};

#endif // CPBRT_CORE_INTEGRATOR_H