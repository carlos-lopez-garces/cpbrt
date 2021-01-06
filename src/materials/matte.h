#include <memory>

#include "core/interaction.h"
#include "core/material.h"
#include "core/memory.h"
#include "core/spectrum.h"

class MatteMaterial : public Material {
private:
    // Texture of diffuse reflectances, also known as diffuse colors or albedos.
    std::shared_ptr<Texture<Spectrum>> Kd;

    // Roughness. When 0, the surface has perfect diffuse Lambertian reflection.
    // When != 0, the surface is microfaceted and its value is the sigma parameter
    // of the Oren-Nayar reflection model (the standard deviation of the microfacet
    // orientation angle).
    std::shared_ptr<Texture<Float>> sigma;

    std::shared_ptr<Texture<Float>> bumpMap;

public:
    MatteMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Float>> &sigma,
        const std::shared_ptr<Texture<Float>> &bumpMap
    ) : Kd(Kd), 
        sigma(sigma), 
        bumpMap(bumpMap)
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