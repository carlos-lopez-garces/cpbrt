#ifndef CPBRT_MATERIALS_DISNEY_H
#define CPBRT_MATERIALS_DISNEY_H

#include "core/cpbrt.h"
#include "core/material.h"
#include "core/paramset.h"

class DisneyMaterial : public Material {
private:
    std::shared_ptr<Texture<Spectrum>> color;

public:
    DisneyMaterial(
        const std::shared_ptr<Texture<Spectrum>> &color
    ) : color(color) {}

    // Evaluates BSDFs at intersection point. PathIntegrator and VolPathIntegrator
    // pass allowMultipleLobes=true. DirectLightingIntegrator passes allowMultipleLobes=false. 
    void ComputeScatteringFunctions(
        SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_MATERIALS_DISNEY_H