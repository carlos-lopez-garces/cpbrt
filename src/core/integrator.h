#ifndef CPBRT_CORE_INTEGRATOR_H
#define CPBRT_CORE_INTEGRATOR_H

#include <memory>

#include "cpbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "film.h"

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

    Film *uiFilm;

protected:
    std::shared_ptr<const Camera> camera;

public:
    SamplerIntegrator(
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler,
        Film *film
    ) : camera(camera), sampler(sampler), uiFilm(film) {}

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

    Spectrum SpecularReflect(
        const RayDifferential &ray,
        const SurfaceInteraction &si,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        int depth
    ) const;
};

#endif // CPBRT_CORE_INTEGRATOR_H