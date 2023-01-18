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
    std::shared_ptr<Texture<Float>> nu;
    std::shared_ptr<Texture<Float>> nv;

    bool remapRoughness;

    // TODO: implement bump mapping.

public:
    SubstrateMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Spectrum>> &Ks,
        const std::shared_ptr<Texture<Float>> &nu,
        const std::shared_ptr<Texture<Float>> &nv,
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness
    )
      : Kd(Kd),
        Ks(Ks),
        nu(nu),
        nv(nv),
        remapRoughness(remapRoughness)
    {}
};

#endif // CPBRT_MATERIALS_SUBSTRATE_H