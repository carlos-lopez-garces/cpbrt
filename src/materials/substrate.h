#ifndef CPBRT_MATERIALS_SUBSTRATE_H
#define CPBRT_MATERIALS_SUBSTRATE_H

#include "core/cpbrt.h"
#include "core/material.h"

// A material with a diffuse underlying surface (diffuse substrate) and a glossy
// specular surface above it.
class SubstrateMaterial : public Material {
private:
    // Diffuse reflectance of the substrate.
    std::shared_ptr<Texture<Spectrum>> Kd;

    // Specular reflectance of the coat.
    std::shared_ptr<Texture<Spectrum>> Ks;

    // Shading normals.
    std::shared_ptr<Texture<Float>> uRoughness;
    std::shared_ptr<Texture<Float>> vRoughness;

    bool remapRoughness;

    // TODO: implement bump mapping.

public:
    SubstrateMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Spectrum>> &Ks,
        const std::shared_ptr<Texture<Float>> &uRough,
        const std::shared_ptr<Texture<Float>> &vRough,
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness
    )
      : Kd(Kd),
        Ks(Ks),
        uRoughness(uRough),
        vRoughness(vRough),
        remapRoughness(remapRoughness)
    {}

    // Evaluates BSDFs at intersection point.
    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_MATERIALS_SUBSTRATE_H