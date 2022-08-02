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