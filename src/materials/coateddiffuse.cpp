#include "coateddiffuse.h"
#include "core/microfacet.h"
#include "core/paramset.h"

CoatedDiffuseMaterial *CreateCoatedDiffuseMaterial(const TextureParams &mp) {
    Spectrum eta = mp.FindSpectrum("eta", 1.5f);
    std::shared_ptr<Texture<Spectrum>> reflectance = mp.GetSpectrumTexture("reflectance", Spectrum(0.5f));
    std::shared_ptr<Texture<Spectrum>> albedo = mp.GetSpectrumTexture("albedo", Spectrum(0.f));
    std::shared_ptr<Texture<Float>> uRoughness = mp.GetFloatTextureOrNull("uroughness");
    std::shared_ptr<Texture<Float>> vRoughness = mp.GetFloatTextureOrNull("vroughness");
    std::shared_ptr<Texture<Float>> thickness = mp.GetFloatTexture("thickness", 0.01f);
    std::shared_ptr<Texture<Float>> g = mp.GetFloatTexture("g", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    int maxDepth = mp.FindInt("maxdepth", 10);
    int nSamples = mp.FindInt("nsamples", 1);
    return new CoatedDiffuseMaterial(reflectance, uRoughness, vRoughness, thickness, albedo, g, eta, bumpMap, remapRoughness, maxDepth, nSamples);
}