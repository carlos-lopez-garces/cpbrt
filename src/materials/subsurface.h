#ifndef CPBRT_MATERIALS_SUBSURFACE_H
#define CPBRT_MATERIALS_SUBSURFACE_H

#include "core/cpbrt.h"
#include "core/material.h"

class SubsurfaceMaterial : public Material {
private:
    // Scales the absorption and scattering coefficients to facilitate change of units.
    const Float scale;

    // Reflectance.
    std::shared_ptr<Texture<Spectrum>> Kr;

    // (Beam) Transmittance.
    std::shared_ptr<Texture<Spectrum>> Kt;

    std::shared_ptr<Texture<Spectrum>> sigma_a;
    std::shared_ptr<Texture<Spectrum>> sigma_s;

    // AlphaX parameter of the Trowbridge-Reitz microfacet distribution (TrowbridgeReitzDistribution).
    std::shared_ptr<Texture<Float>> uRoughness;
    // AlphaY parameter of the Trowbridge-Reitz microfacet distribution (TrowbridgeReitzDistribution).
    std::shared_ptr<Texture<Float>> vRoughness;
    const bool remapRoughness;

    std::shared_ptr<Texture<Float>> bumpMap;

    // Index of refraction of the transmission medium.
    const Float eta;
    
    BSSRDFTable table;

public:
    SubsurfaceMaterial(
        Float scale,
        const std::shared_ptr<Texture<Spectrum>> &Kr,
        const std::shared_ptr<Texture<Spectrum>> &Kt,
        const std::shared_ptr<Texture<Spectrum>> &sigma_a,
        const std::shared_ptr<Texture<Spectrum>> &sigma_s,
        Float g, Float eta,
        const std::shared_ptr<Texture<Float>> &uRoughness,
        const std::shared_ptr<Texture<Float>> &vRoughness,
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness
    ) : scale(scale),
        Kr(Kr),
        Kt(Kt),
        sigma_a(sigma_a),
        sigma_s(sigma_s),
        uRoughness(uRoughness),
        vRoughness(vRoughness),
        bumpMap(bumpMap),
        eta(eta),
        remapRoughness(remapRoughness),
        table(100, 64)
    {
        // TODO: implement.
        ComputeBeamDiffusionBSSRDF(g, eta, &table);
    }

    // Adds BxDFs and BSSRDF to intersection point si.
    void ComputeScatteringFunctions(
        SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_MATERIALS_SUBSURFACE_H