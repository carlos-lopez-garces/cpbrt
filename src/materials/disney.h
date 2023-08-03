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
    // and specular responses; in the diffuse response, it influences diffuse retro-
    // reflection, specifically.
    std::shared_ptr<Texture<Float>> roughness;
    // A soft luster visible at grazing angles, especially on cloth. In [0, 1].
    std::shared_ptr<Texture<Float>> sheen;
    // Controls how much the sheen term's color is tinted by the base color. In [0, 1].
    std::shared_ptr<Texture<Float>> sheenTint;

public:
    DisneyMaterial(
        const std::shared_ptr<Texture<Spectrum>> &color,
        const std::shared_ptr<Texture<Float>> &metallic,
        const std::shared_ptr<Texture<Float>> &roughness,
        const std::shared_ptr<Texture<Float>> &sheen,
        const std::shared_ptr<Texture<Float>> &sheenTint
    ) : color(color),
        metallic(metallic),
        roughness(roughness),
        sheen(sheen),
        sheenTint(sheenTint)
    {}

    // Evaluates BSDFs at intersection point. PathIntegrator and VolPathIntegrator
    // pass allowMultipleLobes=true. DirectLightingIntegrator passes allowMultipleLobes=false. 
    void ComputeScatteringFunctions(
        SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
    ) const;
};

DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_DISNEY_H