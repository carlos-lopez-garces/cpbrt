#include "subsurface.h"
#include "core/microfacet.h"
#include "core/reflection.h"

void SubsurfaceMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
) const {
    // TODO: implement bump mapping.

    // Initialize BSDF.

    Spectrum R = Kr->Evaluate(*si).Clamp();
    Spectrum T = Kt->Evaluate(*si).Clamp();

    // Trowbridge-Reitz microfacet distribution parameters.
    Float uRough = uRoughness->Evaluate(*si);
    Float vRough = vRoughness->Evaluate(*si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    if (R.IsBlack() && T.IsBlack()) {
        // Either one may be black, but not both.
        return;
    }

    bool isSpecular = uRough == 0 && vRough == 0;
    if (isSpecular && allowMultipleLobes) {
        // Not every integrator supports BxDFs that model multiple types of scattering (e.g. those
        // that model both reflection and transmission simultaneously); those who do, pass
        // allowMultipleLobes=true. FresnelSpecularReflectionTransmission models both. If 
        // allowMultipleLobes=false, then we add separately the SpecularReflection and SpecularTransmission
        // BxDFs to the SurfaceInteraction's list of BSDFs.
        si->bsdf->Add(ARENA_ALLOC(arena, FresnelSpecularReflectionTransmission)(R, T, 1.f, eta, mode));
    } else {
        if (remapRoughness) {
            // The alphaX and alphaY parameters of the Trowbridge-Reitz distribution.
            uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
            vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
        }

        MicrofacetDistribution *distribution = nullptr;
        if (!isSpecular) {
          distribution = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(uRough, vRough);
        }
        // Reflection.
        if (!R.IsBlack()) {
            // Evaluates the Fresnel reflectance equation for dielectrics. The index of refraction of the
            // incident medium is 1.0 (i.e. vacuum); the index of refraction of the transmission medium
            // is eta.
            Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);

            if (isSpecular) {
                si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(R, fresnel));
            } else {
                si->bsdf->Add(ARENA_ALLOC(arena, TorranceSparrowMicrofacetReflection)(R, distribution, fresnel));
            }
        }
        // Transmission.
        if (!T.IsBlack()) {
            if (isSpecular) {
                // Incoming medium's index of refraction (1.0, i.e. vacuum) followed by transmission
                // medium's index of refraction.
                si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(T, 1.f, eta, mode));
            } else {
                si->bsdf->Add(ARENA_ALLOC(arena, TorranceSparrowMicrofacetTransmission)(T, distribution, 1.f, eta, mode));
            }
        }
    }

    Spectrum sig_a = scale * sigma_a->Evaluate(*si).Clamp();
    Spectrum sig_s = scale * sigma_s->Evaluate(*si).Clamp();
    si->bssrdf = ARENA_ALLOC(arena, TabulatedBSSRDF)(*si, this, mode, eta, sig_a, sig_s, table);
}

SubsurfaceMaterial *CreateSubsurfaceMaterial(const TextureParams &mp) {
    Float sig_a_rgb[3] = { 0.0011f, 0.0024f, 0.014f };
    Float sig_s_rgb[3] = { 2.55f, 3.21f, 3.77f };

    // Absorption cross-section and scattering coefficient.
    Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb);
    Spectrum sig_s = Spectrum::FromRGB(sig_s_rgb);

    std::string name = mp.FindString("name");
    // See if we have a table of samples of this medium.
    bool found = GetMediumScatteringProperties(name, &sig_a, &sig_s);
    Float g = mp.FindFloat("g", 0.0f);
    if (name != "") {
        if (!found) {
            Warning("Named material \"%s\" not found.  Using defaults.", name.c_str());
        } else {
            // TODO: explain.
            g = 0;
        }
    }

    Float scale = mp.FindFloat("scale", 1.f);
    Float eta = mp.FindFloat("eta", 1.33f);

    std::shared_ptr<Texture<Spectrum>> sigma_a, sigma_s;
    sigma_a = mp.GetSpectrumTexture("sigma_a", sig_a);
    sigma_s = mp.GetSpectrumTexture("sigma_s", sig_s);

    std::shared_ptr<Texture<Spectrum>> Kr = mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    std::shared_ptr<Texture<Spectrum>> Kt = mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    std::shared_ptr<Texture<Float>> uRough = mp.GetFloatTexture("uroughness", 0.f);
    std::shared_ptr<Texture<Float>> vRough = mp.GetFloatTexture("vroughness", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new SubsurfaceMaterial(scale, Kr, Kt, sigma_a, sigma_s, g, eta, uRough, vRough, bumpMap, remapRoughness);
}