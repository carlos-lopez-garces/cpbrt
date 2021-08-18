#include "mirror.h"
#include "core/spectrum.h"
#include "core/reflection.h"
#include "core/paramset.h"
#include "core/texture.h"
#include "core/interaction.h"

void MirrorMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si,
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    // TODO: bump mapping.

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
    Spectrum R = Kr->Evaluate(*si).Clamp();
    if (!R.IsBlack()) {
        // The FresnelNoOp reflectance reflects all incident light in its entirety: no absorption,
        // no transmission; like mirrors do.
        si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(R, ARENA_ALLOC(arena, FresnelNoOp)()));
    }
}