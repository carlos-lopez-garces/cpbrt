#ifndef CPBRT_MATERIALS_PLASTIC_H
#define CPBRT_MATERIALS_PLASTIC_H

#include <memory>

#include "core/cpbrt.h"
#include "core/material.h"

class PlasticMaterial : public Material {
private:
    // Diffuse component. Texture of diffuse reflectances, also known as diffuse
    // colors or albedos.
    std::shared_ptr<Texture<Spectrum>> Kd;

    // Glossy specular component.
    std::shared_ptr<Texture<Spectrum>> Ks;

    // Root mean square slope of microfacets in radians from 0 to Pi/2, measured
    // with respect to the surface normal; main parameter of the microfacet
    // distribution that describes the geometry of the surface.
    //
    // If remapRoughness=true, though, roughness ranges from 0 (perfect specular)
    // to 1 (very rough) and the implementation converts it to root mean square slope.
    std::shared_ptr<Texture<Float>> roughness;
    const bool remapRoughness;

    std::shared_ptr<Texture<Float>> bumpMap;

public:
    PlasticMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Spectrum>> &Ks,
        const std::shared_ptr<Texture<Float>> &roughness,
        const std::shared_ptr<Texture<Float>> &bumpMap,
        bool remapRoughness 
    ) : Kd(Kd),
        Ks(Ks),
        roughness(roughness),
        bumpMap(bumpMap),
        remapRoughness(remapRoughness)
    {}

    void ComputeScatteringFunctions(
        // Differential geometry of surface-ray intersection point.
        SurfaceInteraction *si,
        // For allocating BSDFs.
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

PlasticMaterial *CreatePlasticMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_PLASTIC_H