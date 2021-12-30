#include "scene.h"
#include "spectrum.h"

bool Scene::Intersect(const Ray &ray, SurfaceInteraction *si) const {
    return aggregate->Intersect(ray, si);
}

bool Scene::IntersectP(const Ray &ray) const {
    return aggregate->IntersectP(ray);
}

 bool Scene::IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *si, Spectrum *transmittance) const {
    *transmittance = Spectrum(1.f);
    while (true) {
        bool hitSurface = Intersect(ray, si);

        // Accumulate beam transmittance for ray segment.
        if (ray.medium) {
            // Beam transmittance is multiplicative along points on a ray: Tr(p->p'') = Tr(p->p')Tr(p'->p'').
            // This new factor corresponds to the transmittance between the ray's origin and
            // the point of intersection (if any), which corresponds to the ray's tMax, as set
            // by Intersect().
            *transmittance *= ray.medium->Tr(ray, sampler);
        }

        // Initialize next ray segment or terminate transmittance computation.
        if (!hitSurface) {
            return false;
        }
        if (si->primitive->GetMaterial() != nullptr) {
            return true;
        }

        // The intersected primitive doesn't have a surface, it's just the boundary of a
        // region filled with a participating medium (as indicated by the MediumInterface
        // of this optically inactive region). These boundaries don't deflect the ray
        // or change the radiance that they carry in any way.
        //
        // It is assumed that the indices of refraction of the region at the ray's origin
        // and the region the ray is entering are equal, despite being different media.
        // Otherwise, the region would've had a material with a BTDF.
        ray = si->SpawnRay(ray.d);
    }
 }