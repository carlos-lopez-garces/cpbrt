#ifndef CPBRT_LIGHTS_DISTANT_H
#define CPBRT_LIGHTS_DISTANT_H

#include "core/cpbrt.h"
#include "core/light.h"
#include "core/medium.h"

class DistantLight : public Light {
private:
    const Spectrum L;
    const Vector3f wLight;

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
};

#endif // CPBRT_LIGHTS_DISTANT_H