#include "distant.h"
#include "core/paramset.h"

Spectrum DistantLight::Sample_Li(
    const Interaction &it,
    const Point2f &u,
    Vector3f *wi,
    Float *pdf,
    VisibilityTester *vis
) const {
    // Constant incident direction.
    *wi = wLight;
    *pdf = 1;

    // Place p1 outside the scene along the light source's direction. A distant
    // light doesn't emit radiance from any particular location, just along the
    // same direction.
    Point3f pBeyondScene = it.p + wLight * (2 * worldRadius);
    *vis = VisibilityTester(it, Interaction(pBeyondScene, it.time, mediumInterface));

    // Constant incident radiance.
    return L;
}

// Power, aka radiant flux, is the total amount of energy passing through a surface
// per unit time. For distant lights, power is Phi = AL, where L is emitted radiance
// and A is total *unoccluded* surface area.
//
// It's impractical to compute A exactly, so it is approximated by the area of a disk
// inside the scene's bounding sphere with a normal vector parallel to the direction
// of the light.
Spectrum DistantLight::Power() const {
    return L * Pi * worldRadius * worldRadius;
}

std::shared_ptr<DistantLight> CreateDistantLight(
    const Transform &light2world,
    const ParamSet &paramSet
) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum scale = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point3f from = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
    Point3f to = paramSet.FindOnePoint3f("to", Point3f(0, 0, 1));
    Vector3f direction = from - to;
    return std::make_shared<DistantLight>(light2world, L*scale, direction);
}