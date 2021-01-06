#include "core/reflection.h"
#include "plastic.h"

void PlasticMaterial::ComputeScatteringFunctions(
    // Differential geometry of surface-ray intersection point.
    SurfaceInteraction *si,
    // For allocating BSDFs.
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    // TODO: Perform bump mapping with bumpMap, if present, to compute shading normal.

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Initialize diffuse component. Reflectance aka diffuse color aka albedo.
    Spectrum kd = Kd->Evaluate(*si).Clamp();
    if (!kd.IsBlack()) {
        si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));
    }

    // Initialize glossy specular component.
    Spectrum ks = Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack()) {
        // TODO: implement when glossy specular BxDFs are ready.
    }
}