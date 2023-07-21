#include "disney.h"
#include "core/interaction.h"
#include "core/reflection.h"

void DisneyMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
) const {
    // TODO: bump mapping.

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Compute diffuse component.
    
    Spectrum colour = color->Evaluate(*si).Clamp();
    Float metal = metallic->Evaluate(*si);
    Float IOR = eta->Evaluate(*si);

    si->bsdf->Add(ARENA_ALLOC(arena, DisneyDiffuseReflection)(
        colour
    ));

    si->bsdf->Add(ARENA_ALLOC(arena, LambertianTransmission)(
        colour
    ));
}