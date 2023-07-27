#ifndef CPBRT_MATERIALS_DISNEY_H
#define CPBRT_MATERIALS_DISNEY_H

#include "core/cpbrt.h"
#include "core/material.h"
#include "core/paramset.h"

class DisneyMaterial : public Material {
private:
    std::shared_ptr<Texture<Spectrum>> color;
    // In [0, 1], where 0 = dielectric and 1 = metallic. (Most parameters are in [0,1].)
    std::shared_ptr<Texture<Float>> metallic;
    // Degree of roughness, in [0, 1], where 0 is smooth. Influences both the diffuse
    // and specular responses.
    std::shared_ptr<Texture<Float>> roughness;
    // Index of refraction IOR of the inside of the object.
    std::shared_ptr<Texture<Float>> eta;

public:
    DisneyMaterial(
        const std::shared_ptr<Texture<Spectrum>> &color,
        const std::shared_ptr<Texture<Float>> &metallic,
        const std::shared_ptr<Texture<Float>> &roughness
    ) : color(color),
        metallic(metallic),
        roughness(roughness)
    {}

    // Evaluates BSDFs at intersection point. PathIntegrator and VolPathIntegrator
    // pass allowMultipleLobes=true. DirectLightingIntegrator passes allowMultipleLobes=false. 
    void ComputeScatteringFunctions(
        SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
    ) const;
};

DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_DISNEY_H