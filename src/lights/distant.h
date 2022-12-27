#ifndef CPBRT_LIGHTS_DISTANT_H
#define CPBRT_LIGHTS_DISTANT_H

#include "core/cpbrt.h"
#include "core/light.h"
#include "core/medium.h"
#include "core/scene.h"

class DistantLight : public Light {
private:
    // Emitted radiance.
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

    // This kind of light source is not sampled: the incident direction w_i is the light
    // source's direction with probability 1. The magnitude of L_i is constant too.
    Spectrum Sample_Li(
        const Interaction &it,
        // A random variable is not needed because there's no random sampling done here.
        // The incident direction and radiance are constant.
        const Point2f &u,
        // Always wLight.
        Vector3f *wi,
        // Always 1.
        Float *pdf,
        // Visibility testing is done by placing a point outside the world/scene bounds
        // along the light source's constant direction.
        VisibilityTester *vis
    ) const override;

    // Evaluates the PDF of the distribution of directions of this light source for
    // the input incident direction. The distribution in this case is constant: wLight
    // with probability 1.
    Float Pdf_Li(const Interaction &it, const Vector3f &wi) const;

    // Power, aka radiant flux, is the total amount of energy passing through a surface
    // per unit time. For distant lights, power is Phi = AL, where L is emitted radiance
    // and A is total *unoccluded* surface area.  
    Spectrum Power() const override;
};

std::shared_ptr<DistantLight> CreateDistantLight(
    const Transform &light2world,
    const ParamSet &paramSet
);

#endif // CPBRT_LIGHTS_DISTANT_H