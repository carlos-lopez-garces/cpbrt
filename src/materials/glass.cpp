#include "glass.h"
#include "core/error.h"
#include "core/interaction.h"
#include "core/paramset.h"
#include "core/reflection.h"
#include "core/spectrum.h"
#include "core/texture.h"

void GlassMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
) const {
    // IOR of the other side of the boundary.
    Float eta = etas->Evaluate(*si);

    // TODO: evaluate uRoughness and vRoughness when MicrofacetTransmission is implemented.
    Float uRough = 0.f;
    Float vRough = 0.f;

    // Reflectance.
    Spectrum R = Kr->Evaluate(*si).Clamp();
    // Transmittance.
    Spectrum T = Kt->Evaluate(*si).Clamp();

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    if (R.IsBlack() && T.IsBlack()) {
        return;   
    }

    bool isSpecular = uRough == 0. && vRough == 0.;
    if (isSpecular && allowMultipleLobes) {
        // PathIntegrator and VolPathIntegrator pass allowMultipleLobes=true.
        si->bsdf->Add(ARENA_ALLOC(arena, FresnelSpecularReflectionTransmission)(R, T, 1., eta, mode));
    } else {
        // DirectLightingIntegrator passes allowMultipleLobes=false.

        if (!R.IsBlack()) {
            Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
            if (isSpecular) {
                si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(R, fresnel));
            } else {
                // TODO: add MicrofacetReflection after implementing MicrofacetTransmission.
                Error("Glass with MicrofacetReflection is not implemented yet.");
            }
        }

        if (!T.IsBlack()) {
            if (isSpecular) {
                si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(T, 1.f, eta, mode));
            } else {
                // TODO: implement MicrofacetTransmission.
                Error("Glass with MicrofacetTransmission is not implemented yet.");
            } 
        }
    }
}

GlassMaterial *CreateGlassMaterial(const TextureParams &mp) {
    // Reflectance.
    std::shared_ptr<Texture<Spectrum>> Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    // Transmittance.
    std::shared_ptr<Texture<Spectrum>> Kt = mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    // Index of refraction IOR of the inside of the object.
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTextureOrNull("eta");
    if (!eta) eta = mp.GetFloatTexture("index", 1.5f);
    std::shared_ptr<Texture<Float>> roughu = mp.GetFloatTexture("uroughness", 0.f);
    std::shared_ptr<Texture<Float>> roughv = mp.GetFloatTexture("vroughness", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new GlassMaterial(Kr, Kt, roughu, roughv, eta, bumpMap, remapRoughness);
}