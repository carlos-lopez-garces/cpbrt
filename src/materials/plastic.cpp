#include "plastic.h"
#include "core/spectrum.h"
#include "core/reflection.h"
#include "core/paramset.h"
#include "core/texture.h"
#include "core/interaction.h"
#include "core/microfacet.h"

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
        // etaI = 1.0 and etaT = 1.5 are the indices of refraction for plastic.
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, 1.5f);

        // Create Beckmann-Spizzichino microfacet distribution. Can also use Trowbridge-Reitz,
        // but it isn't implemented yet.
        Float rough = roughness->Evaluate(*si);
        if (remapRoughness) {
            // Roughness was specified in the range [0,1], instead of as the root mean square
            // (RMS) slope of the microfacets, alpha. Convert to RMS slope.
            rough = BeckmannDistribution::RoughnessToAlpha(rough);
        }
        MicrofacetDistribution *distribution = ARENA_ALLOC(arena, BeckmannDistribution)(rough, rough);

        BxDF *specularBxDF = ARENA_ALLOC(arena, TorranceSparrowMicrofacetReflection)(ks, distribution, fresnel);
        si->bsdf->Add(specularBxDF);
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