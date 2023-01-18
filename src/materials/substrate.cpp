#include "substrate.h"

void SubstrateMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si,
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    // TODO: bump mapping.

    // Reflectances of diffuse substrate and glossy coat.
    Spectrum kd = Kd->Evaluate(*si).Clamp();
    Spectrum ks = Ks->Evaluate(*si).Clamp();

    // Roughness.
    Float uRough = uRoughness->Evaluate(*si);
    Float vRough = vRoughness->Evaluate(*si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    if (!d.IsBlack() || !s.IsBlack()) {
        if (remapRoughness) {
            uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
            vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
        }

        MicrofacetDistribution *distribution = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(uRough, vRough);

        si->bsdf->Add(ARENA_ALLOC(arena, AshikhminShirleyReflection)(kd, ks, distrib));
    }
}