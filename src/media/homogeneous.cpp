#include "homogeneous.h"
#include "core/memory.h"

Spectrum HomogeneousMedium::Tr(const Ray &ray, Sampler &sampler) const {
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}

Spectrum HomogeneousMedium::Sample(
    const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi
) const {
    // Sample a channel and distance along the ray. 
    int channel = std::min((int)(sampler.Get1D() * Spectrum::nSamples), Spectrum::nSamples - 1);
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t[channel];
    Float t = std::min(dist * ray.d.Length(), ray.tMax);
    bool sampledMedium = t < ray.tMax;
    if (sampledMedium) {
        *mi = MediumInteraction(ray(t), -ray.d, ray.time, this, ARENA_ALLOC(arena, HenyeyGreensteinPhaseFunction)(g));
    }

    // Compute the transmittance and sampling density. 
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());

    // Return weighting factor for scattering from homogeneous medium. 
    Spectrum density = sampledMedium ? (sigma_t * Tr) : Tr;
    Float pdf = 0;
    for (int i = 0; i < Spectrum::nSamples; ++i) {
        pdf += density[i];
    }
    pdf *= 1 / (Float)Spectrum::nSamples;
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}