#ifndef CPBRT_MATERIALS_GLASS_H
#define CPBRT_MATERIALS_GLASS_H

#include "core/cpbrt.h"
#include "core/material.h"

class GlassMaterial : public Material {
private:
    // Reflectance.
    std::shared_ptr<Texture<Spectrum>> Kr;
    // Transmittance.
    std::shared_ptr<Texture<Spectrum>> Kt;
    // Adjusts the microfacet distribution's roughness.
    // TODO: implement MicrofacetTransmission.
    std::shared_ptr<Texture<Float>> uRoughness;
    std::shared_ptr<Texture<Float>> vRoughness;
    // Index of refraction IOR of the inside of the object.
    std::shared_ptr<Texture<Float>> etas;
    // TODO: implement bump mapping.
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;

public:
    GlassMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kr,
        const std::shared_ptr<Texture<Spectrum>> &Kt,
        const std::shared_ptr<Texture<Float>> &uRoughness,
        const std::shared_ptr<Texture<Float>> &vRoughness,
        const std::shared_ptr<Texture<Float>> &index,
        // TODO: implement bump mapping.
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness
    ) : Kr(Kr),
        Kt(Kt),
        uRoughness(uRoughness),
        vRoughness(vRoughness),
        index(index),
        bumpMap(bumpMap),
        remapRoughness(remapRoughness)
    {}

    // Evaluates BSDFs at intersection point. PathIntegrator and VolPathIntegrator
    // pass allowMultipleLobes=true. DirectLightingIntegrator passes allowMultipleLobes=false. 
    void ComputeScatteringFunctions(
        SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_MATERIALS_GLASS_H