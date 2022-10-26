#include "ward.h"
#include "core/paramset.h"
#include "core/reflection.h"
#include "core/interaction.h"
#include "core/texture.h"
#include "core/interaction.h"

void WardMaterial::ComputeScatteringFunctions(
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
    Spectrum ks = Ks->Evaluate(*si).Clamp();

    Float alphax = Clamp(alphaX->Evaluate(*si), 0, 1);
    Float alphay = Clamp(alphaY->Evaluate(*si), 0, 1);

    if (!kd.IsBlack()) {
        si->bsdf->Add(ARENA_ALLOC(arena, WardReflection)(kd, ks, alphax, alphay));
    }
}

WardMaterial *CreateWardMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    std::shared_ptr<Texture<Spectrum>> Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> alphaX = mp.GetFloatTexture("alphaX", 1.f);
    std::shared_ptr<Texture<Float>> alphaY = mp.GetFloatTexture("alphaY", 1.f);

    return new WardMaterial(Kd, Ks, alphaX, alphaY);
}