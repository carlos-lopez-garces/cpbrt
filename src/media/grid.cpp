#include "grid.h"
#include "core/cpbrt.h"
#include "core/geometry.h"

Float GridDensityMedium::Density(const Point3f &p) const {
    // Compute voxel coordinate and offsets for p.
    Point3f pSamples(p.x*nx - 0.5f, p.y*ny - 0.5f, p.z*nz - 0.5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // Trilinearly interpolate density values to compute local density.
    Float d00 = Lerp(d.x, D(pi), D(pi+Vector3i(1,0,0)));
    Float d10 = Lerp(d.x, D(pi+Vector3i(0,1,0)), D(pi+Vector3i(1,1,0)));
    Float d01 = Lerp(d.x, D(pi+Vector3i(0,0,1)), D(pi+Vector3i(1,0,1)));
    Float d11 = Lerp(d.x, D(pi+Vector3i(0,1,1)), D(pi+Vector3i(1,1,1)));
    Float d0 = Lerp(d.y, d00, d10);
    Float d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

Spectrum GridDensityMedium::Sample(
    const Ray &rWorld, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi
) const {
    Ray ray = WorldToMedium(Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));

    // Compute [t_min, t_max] interval of ray's overlap with medium bounds.
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) {
        return Spectrum(1.f);
    }

    // Run delta-tracking iterations to sample a medium interaction.
    Float t = tMin;
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * reciprocalMaxDensity / sigma_t;  
        if (t >= tMax) {
            break;
        }
        if (Density(ray(t)) * reciprocalMaxDensity > sampler.Get1D()) {
            // Populate mi with medium interaction information and return.
            PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreensteinPhaseFunction)(g);
            *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this, phase);
            return sigma_s / sigma_t;
        }
    }
    return Spectrum(1.f);
}