#include "substrate.h"
#include "core/interaction.h"
#include "core/reflection.h"
#include "core/paramset.h"

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

    if (!kd.IsBlack() || !ks.IsBlack()) {
        if (remapRoughness) {
            uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
            vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
        }

        MicrofacetDistribution *distribution = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(uRough, vRough);

        si->bsdf->Add(ARENA_ALLOC(arena, AshikhminShirleyReflection)(kd, ks, distribution));
    }
}

SubstrateMaterial *CreateSubstrateMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(.5f));
    std::shared_ptr<Texture<Spectrum>> Ks = mp.GetSpectrumTexture("Ks", Spectrum(.5f));
    std::shared_ptr<Texture<Float>> uRoughness = mp.GetFloatTexture("uroughness", .1f);
    std::shared_ptr<Texture<Float>> vRoughness = mp.GetFloatTexture("vroughness", .1f);
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new SubstrateMaterial(Kd, Ks, uRoughness, vRoughness, remapRoughness);
}
