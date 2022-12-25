#ifndef CPBRT_LIGHTS_DISTANT_H
#define CPBRT_LIGHTS_DISTANT_H

#include "core/cpbrt.h"
#include "core/light.h"
#include "core/medium.h"

class DistantLight : public Light {
private:
    const Spectrum L;

    // Direction in world space. Radiance emitted by this light source arrives at
    // surfaces in parallel beams in this direction. 
    const Vector3f wLight;

    Point3f worldCenter;
    Float worldRadius;

public:
    // Unlike PointLight and DiffuseAreaLight, DistantLight cannot be embedded in
    // participating media because it is supposed to be very far away. The absorption
    // effects that any participating media at the source location might have is
    // assumed to have taken place already.
    DistantLight(
        const Transform &LightToWorld,
        const Spectrum &L,
        const Vector3f &wLight
    ) : Light((int) LightFlags::DeltaDirection, LightToWorld, MediumInterface()),
        L(L),
        wLight(Normalize(LightToWorld(wLight)))
    {}

    // Get world bounds.
    void Preprocess(const Scene &scene) {
        scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    }
};

#endif // CPBRT_LIGHTS_DISTANT_H