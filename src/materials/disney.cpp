#include "disney.h"
#include "core/interaction.h"
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
    // Float IOR = eta->Evaluate(*si);

    si->bsdf->Add(ARENA_ALLOC(arena, DisneyDiffuseReflection)(
        colour * diffuseWeight
    ));
}

DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> color = mp.GetSpectrumTexture("color", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> metallic = mp.GetFloatTexture("metallic", 0.f);
    // std::shared_ptr<Texture<Float>> eta = mp.GetFloatTexture("eta", 1.5f);
    return new DisneyMaterial(color, metallic);
}