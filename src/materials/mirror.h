#ifndef CPBRT_MATERIALS_MIRROR_H
#define CPBRT_MATERIALS_MIRROR_H

#include "core/cpbrt.h"
#include "core/material.h"

// The BRDF of a mirror is one that reflects all the incident light in its
// entirety, as opposed to just a fraction of it like other BRDFs do.
class MirrorMaterial : public Material {
private:
    std::shared_ptr<Texture<Spectrum>> Kr;
    
    // TODO: add bump map.

public:
    MirrorMaterial(const std::shared_ptr<Texture<Spectrum>> &r) {
        Kr = r;
    }

    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

MirrorMaterial *CreateMirrorMaterial(const TextureParams &mp);

#endif // CPBRT_MATERIALS_MIRROR_H