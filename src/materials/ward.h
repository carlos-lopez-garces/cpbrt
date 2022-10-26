#ifndef CPBRT_MATERIALS_WARD_H
#define CPBRT_MATERIALS_WARD_H

#include <memory>

#include "core/cpbrt.h"
#include "core/material.h"

class WardMaterial : public Material {
private:
    std::shared_ptr<Texture<Spectrum>> Kd;
    std::shared_ptr<Texture<Spectrum>> Ks;
    std::shared_ptr<Texture<Float>> alphaX;
    std::shared_ptr<Texture<Float>> alphaY;


public:
    WardMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kd,
        const std::shared_ptr<Texture<Spectrum>> &Ks,
        const std::shared_ptr<Texture<Float>> &alphaX,
        const std::shared_ptr<Texture<Float>> &alphaY
    ) : Kd(Kd),
        Ks(Ks),
        alphaX(alphaX),
        alphaY(alphaY)
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

WardMaterial *CreateWardMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_WARD_H