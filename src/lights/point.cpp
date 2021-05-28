#include "point.h"
#include "core/scene.h"
#include "core/paramset.h"
#include "core/sampling.h"

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

std::shared_ptr<PointLight> CreatePointLight(
    const Transform &light2world,
    const Medium *medium,
    const ParamSet &paramSet
) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point3f P = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
    Transform l2w = Translate(Vector3f(P.x, P.y, P.z)) * light2world;
    return std::make_shared<PointLight>(l2w, medium, I * sc);
}