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
        si->bsdf->Add(arena, FresnelSpecularReflectionTransmission)(R, T, 1., eta, mode);
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