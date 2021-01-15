#include "core/interaction.h"
#include "point.h"

Spectrum PointLight::Sample_Li(
    const Interaction &it,
    const Point2f &u,
    Vector3f *wi,
    Float *pdf,
    VisibilityTester *vis
) const {
    // A PointLight can only illuminate a given point from a single incident direction.
    *wi = Normalize(pLight - it.p);

    // TODO: explain.
    *pdf = 1.f; 

    *vis = VisibilityTester(it, Interaction(pLight, it.time, mediumInterface));

    // Flux (power) along a given direction falls off with distance at the rate of 1/m^2.
    return I / DistanceSquared(pLight, it.p);
}

Spectrum PointLight::Power() const {
    return 4 * Pi * I;
}