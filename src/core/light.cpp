#include "light.h"
#include "scene.h"
#include "sampling.h"
#include "paramset.h"

Spectrum Light::Le(const RayDifferential &rd) const {
    // This is typically implemented only by environment lights, which are
    // sampled by rays that escape the bounds of the scene without hitting
    // anything.
    return Spectrum(0.f);
}

bool VisibilityTester::Unoccluded(const Scene &scene) const {
    // TODO: implement Interaction::SpawnRayTo.
    return !scene.IntersectP(p0.SpawnRayTo(p1));
}