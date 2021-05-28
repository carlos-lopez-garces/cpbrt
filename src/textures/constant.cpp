#include "textures/constant.h"

ConstantTexture<Float> *CreateConstantFloatTexture(
    const Transform &tex2world,
    const TextureParams &tp
) {
    return new ConstantTexture<Float>(tp.FindFloat("value", 1.f));
}

ConstantTexture<Spectrum> *CreateConstantSpectrumTexture(
    const Transform &tex2world, 
    const TextureParams &tp
) {
    return new ConstantTexture<Spectrum>(tp.FindSpectrum("value", Spectrum(1.f)));
}