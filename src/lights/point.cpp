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

    // A PointLight has a delta distribution of direction; it illuminates a given
    // point of incidence from a single direction with probability 1.
    *pdf = 1.f; 

    *vis = VisibilityTester(it, Interaction(pLight, it.time, mediumInterface));

    // Flux (power) along a given direction falls off with distance at the rate of 1/m^2.
    return I / DistanceSquared(pLight, it.p);
}

Float PointLight::Pdf_Li(const Interaction &it, const Vector3f &wi) const {
    // There's no chance that the given point of incidence will be illuminated by this
    // light source from the given (random) incident direction, unless it is the one
    // computed by Sample_Li.
    return 0.f;
}

Spectrum PointLight::Power() const {
    return 4 * Pi * I;
}