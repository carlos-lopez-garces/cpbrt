#ifndef CPBRT_MATERIALS_COATED_DIFFUSE_H
#define CPBRT_MATERIALS_COATED_DIFFUSE_H

#include "core/cpbrt.h"
#include "core/material.h"
#include "core/spectrum.h"

// CoatedDiffuseMaterial describes scattering from ..., where ....
class CoatedDiffuseMaterial : public Material {
private:
    // Index of refraction.
    Spectrum eta;
    std::shared_ptr<Texture<Spectrum>> reflectance;
    std::shared_ptr<Texture<Spectrum>> albedo;
    // Adjusts the microfacet distribution's roughness.
    std::shared_ptr<Texture<Float>> uRoughness;
    std::shared_ptr<Texture<Float>> vRoughness;
    std::shared_ptr<Texture<Float>> thickness;
    std::shared_ptr<Texture<Float>> g;
    // TODO: implement bump mapping.
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
    int maxDepth;
    int nSamples;

public:
    CoatedDiffuseMaterial(
        const std::shared_ptr<Texture<Spectrum>> &reflectance,
        const std::shared_ptr<Texture<Float>> &uRoughness,
        const std::shared_ptr<Texture<Float>> &vRoughness,
        const std::shared_ptr<Texture<Float>> &thickness,
        const std::shared_ptr<Texture<Spectrum>> &albedo,
        const std::shared_ptr<Texture<Float>> &g,
        const Spectrum &eta,
        // TODO: implement bump mapping.
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness,
        int maxDepth,
        int nSamples
    ) : reflectance(reflectance),
        uRoughness(uRoughness),
        vRoughness(vRoughness),
        thickness(thickness),
        albedo(albedo),
        g(g),
        eta(eta),
        bumpMap(bumpMap),
        remapRoughness(remapRoughness),
        maxDepth(maxDepth),
        nSamples(nSamples)
    {}

    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

CoatedDiffuseMaterial *CreateCoatedDiffuseMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_COATED_DIFFUSE_H