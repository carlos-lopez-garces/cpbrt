#ifndef CPBRT_MATERIALS_METAL_H
#define CPBRT_MATERIALS_METAL_H

#include "core/cpbrt.h"
#include "core/material.h"
#include "core/spectrum.h"

class MetalMaterial : public Material {
private:
    std::shared_ptr<Texture<Spectrum>> eta;
    std::shared_ptr<Texture<Spectrum>> k;
    std::shared_ptr<Texture<Float>> roughness;
    std::shared_ptr<Texture<Float>> uRoughness;
    std::shared_ptr<Texture<Float>> vRoughness;
    // TODO: implement bump mapping.
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;

public:
    MetalMaterial(
        const std::shared_ptr<Texture<Spectrum>> &eta,
        const std::shared_ptr<Texture<Spectrum>> &k,
        const std::shared_ptr<Texture<Float>> &roughness,
        const std::shared_ptr<Texture<Float>> &uRoughness,
        const std::shared_ptr<Texture<Float>> &vRoughness,
        // TODO: implement bump mapping.
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness
    ) : eta(eta), 
        k(k),
        roughness(roughness),
        uRoughness(uRoughness),
        vRoughness(vRoughness),
        bumpMap(bumpMap),
        remapRoughness(remapRoughness)
    {}

    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_MATERIALS_METAL_H