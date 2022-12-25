#include "distant.h"

Spectrum DistantLight::Sample_Li(
    const Interaction &it,
    const Point2f &u,
    Vector3f *wi,
    Float *pdf,
    VisibilityTester *vis
) {
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