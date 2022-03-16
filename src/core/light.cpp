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
    Ray ray(p0.SpawnRayTo(p1));
    Spectrum Tr(1.f);
    while (true) {
        SurfaceInteraction si;
        bool hitSurface = scene.Intersect(ray, &si);
        if (hitSurface && si.primitive->GetMaterial() != nullptr) {
            // The surface is opaque. No transmittance.
            return Spectrum(0.0f);
        }
        if (ray.medium) {
            Tr *= ray.medium->Tr(ray, sampler);
        }
        if (!hitSurface) {
            break;
        }
        ray = si.SpawnRayTo(p1);
    }
    return Tr;
}