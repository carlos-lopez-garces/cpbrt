#include "disney.h"
#include "core/interaction.h"
#include "core/microfacet.h"
#include "core/reflection.h"

void DisneyMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes
) const {
    // TODO: bump mapping.

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Compute diffuse component.
    
    Spectrum colour = color->Evaluate(*si).Clamp();
    Float metallicWeight = metallic->Evaluate(*si);
    Float diffuseWeight = 1.f - metallicWeight;
    Float rough = roughness->Evaluate(*si);
    // The CIE Y curve is luminance.
    Float luminance = colour.y();
    // TODO: explain?
    Spectrum colourTint = luminance > 0 ? colour / luminance : Spectrum(1.f);

    if (diffuseWeight) {
        // DisneyDiffuseReflection implements the diffuse term of the f_d equation (4) in
        // Extending the Disney BRDF to a BSDF with Integrated Subsurface Scattering by Burley.
        // The f_retro-reflection term of f_d is implemented in DisneyRetroReflection.
        si->bsdf->Add(ARENA_ALLOC(arena, DisneyDiffuseReflection)(
            // The diffuse response decreases as the metallic parameter increases.
            colour * diffuseWeight
        ));

        // DisneyRetroReflection implements the f_retro-reflection term of the f_d equation (4) in
        // Extending the Disney BRDF to a BSDF with Integrated Subsurface Scattering by Burley.
        si->bsdf->Add(ARENA_ALLOC(arena, DisneyRetroReflection)(
            colour * diffuseWeight,
            rough
        ));

        // TODO: didn't make a difference.
        Float sheenWeight = sheen->Evaluate(*si);
        if (sheenWeight > 0.f) {
            // The larger the tint weight, the stronger the influence of the base color is on
            // the sheen color.
            Float sheenTintWeight = sheenTint->Evaluate(*si);
            Spectrum sheenColour = Lerp(sheenTintWeight, Spectrum(1.f), colourTint);
            si->bsdf->Add(ARENA_ALLOC(arena, DisneySheenReflection)(
                // sheenColour * sheenWeight * diffuseWeight
                Spectrum(1.f) 
            ));
        }
    }
}

DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> color = mp.GetSpectrumTexture("color", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> metallic = mp.GetFloatTexture("metallic", 0.f);
    std::shared_ptr<Texture<Float>> roughness = mp.GetFloatTexture("roughness", 0.5f);
    std::shared_ptr<Texture<Float>> sheen = mp.GetFloatTexture("sheen", 0.f);
    std::shared_ptr<Texture<Float>> sheenTint = mp.GetFloatTexture("sheenTint", 0.5f);
    return new DisneyMaterial(color, metallic, roughness, sheen, sheenTint);
}