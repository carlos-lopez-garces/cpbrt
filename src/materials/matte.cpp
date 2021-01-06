#include "core/reflection.h"
#include "matte.h"

void MatteMaterial::ComputeScatteringFunctions(
    // Differential geometry of surface-ray intersection point.
    SurfaceInteraction *si,
    // For allocating BSDFs.
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    // TODO: Perform bump mapping with bumpMap, if present, to compute shading normal.

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Reflectance aka diffuse color aka albedo. Clamp to [0, infinity).
    Spectrum kd = Kd->Evaluate(*si).Clamp();

    // Oren-Nayar diffuse reflection model parameter.
    Float sigma = Clamp(sigma->Evaluate(*si), 0, 90);

    if (!kd.IsBlack()) {
        if (sigma == 0) {
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));
        } else {
            // TODO: implement OrenNayarReflection.
            si->bsdf->Add(ARENA_ALLOC(arena, OrenNayarReflection)(kd, sigma));
        }
    }
}