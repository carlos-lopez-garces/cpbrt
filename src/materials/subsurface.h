#ifndef CPBRT_MATERIALS_SUBSURFACE_H
#define CPBRT_MATERIALS_SUBSURFACE_H

#include "core/cpbrt.h"
#include "core/material.h"

class SubsurfaceMaterial : public Material {
private:
    // Scales the absorption and scattering coefficients to facilitate change of units.
    const Float scale;
    std::shared_ptr<Texture<Spectrum>> Kr;
    std::shared_ptr<Texture<Spectrum>> Kt;
    std::shared_ptr<Texture<Spectrum>> sigma_a;
    std::shared_ptr<Texture<Spectrum>> sigma_s;
    std::shared_ptr<Texture<Float>> uRoughness;
    std::shared_ptr<Texture<Float>> vRoughness;
    std::shared_ptr<Texture<Float>> bumpMap;
    const Float eta;
    const bool remapRoughness;
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
        ComputeBeamDiffusionBSSRDF(g, eta, &table);
    }
};

#endif // CPBRT_MATERIALS_SUBSURFACE_H