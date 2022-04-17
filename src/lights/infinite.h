#ifndef CPBRT_LIGHTS_INFINITE_H
#define CPBRT_LIGHTS_INFINITE_H

#include "core/light.h"
#include "core/mipmap.h"
#include "core/scene.h"
#include "core/spectrum.h"
#include "core/transform.h"

class InfiniteAreaLight : public Light {
private:
    std::unique_ptr<MIPMap<RGBSpectrum>> LMap;
    Point3f worldCenter;
    Float worldRadius;

public:
    // If the light is backed by a texture, each texel's value gets multiplied
    // by power; otherwise, the light uses a single texel with power for its value.
    InfiniteAreaLight(
        const Transform &lightToWorld,
        const Spectrum &power,
        int nSamples,
        const std::string &textureFilepath
    );

    void Preprocess(const Scene &scene) {
        scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    }

    // Computes total power emitted.
    Spectrum Power() const;
};

#endif // CPBRT_LIGHTS_INFINITE_H