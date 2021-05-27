#ifndef CPBRT_TEXTURES_CONSTANT_H
#define CPBRT_TEXTURES_CONSTANT_H

#include "core/interaction.h"
#include "core/texture.h"

template <typename T> class ConstantTexture : public Texture<T> {
private:
    // Constant.  
    T value;

public:
    ConstantTexture(const T &value) : value(value) {}

    T Evaluate(const SurfaceInteraction &) const {
        return value;
    }
};

ConstantTexture<Float> *CreateConstantFloatTexture(const Transform &tex2world, const TextureParams &tp);

ConstantTexture<Spectrum> *CreateConstantSpectrumTexture(const Transform &tex2world, const TextureParams &tp);

#endif // CPBRT_TEXTURES_CONSTANT_H