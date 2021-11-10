#include "matte.h"
#include "core/paramset.h"
#include "core/reflection.h"
#include "core/interaction.h"
#include "core/texture.h"
#include "core/interaction.h"

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

    // Oren-Nayar microfaceted diffuse reflection model parameter: the standard deviation
    // of the orientation angle of the microfacets.
    Float sig = Clamp(sigma->Evaluate(*si), 0, 90);

    if (!kd.IsBlack()) {
        if (sig == 0) {
            // Perfect Lambertian reflection. Adequate for surfaces that aren't too rough.
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));
        } else {
            // Oren-Nayar reflection. Models reflection from rough surfaces.
            si->bsdf->Add(ARENA_ALLOC(arena, OrenNayarReflection)(kd, sig));
        }
    }
}

MatteMaterial *CreateMatteMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    return new MatteMaterial(Kd, sigma, bumpMap);
}