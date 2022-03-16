#ifndef CPBRT_MATERIALS_METAL_H
#define CPBRT_MATERIALS_METAL_H

#include "core/cpbrt.h"
#include "core/material.h"
#include "core/spectrum.h"

// MetalMaterial describes scattering from metals, where the index of refraction (eta)
// and the absorption coefficient (k) describe metals' reflectance spectra. These and 
// a roughness parameter, which adjusts the microfacet distributions roughness, describe
// the overall material. 
class MetalMaterial : public Material {
private:
    // Index of refraction.
    std::shared_ptr<Texture<Spectrum>> eta;
    // Absorption coefficient.
    std::shared_ptr<Texture<Spectrum>> k;
    // Adjusts the microfacet distributions roughness.
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

MetalMaterial *CreateMetalMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_METAL_H