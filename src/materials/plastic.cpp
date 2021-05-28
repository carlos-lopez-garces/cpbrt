#include "plastic.h"
#include "core/spectrum.h"
#include "core/reflection.h"
#include "core/paramset.h"
#include "core/texture.h"
#include "core/interaction.h"

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

PlasticMaterial *CreatePlasticMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Float>> roughness = mp.GetFloatTexture("roughness", .1f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new PlasticMaterial(Kd, Ks, roughness, bumpMap, remapRoughness);
}