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
    return !scene.IntersectP(p0.SpawnRayTo(p1));
}

Spectrum VisibilityTester::Tr(const Scene &scene, Sampler &sampler) const {
    // TODO: implement when implementing volume scattering.
    Spectrum Tr(0.f);
    return Tr;
}