#ifndef CPBRT_CORE_SCENE_H
#define CPBRT_CORE_SCENE_H

#include <memory>

#include "cpbrt.h"
#include "geometry.h"
#include "primitive.h"
#include "light.h"

class Scene {
private:
    // All the primitives in the scene.
    std::shared_ptr<Primitive> aggregate;

    Bounds3f worldBound;

public:
    std::vector<std::shared_ptr<Light>> lights;

    Scene(
        std::shared_ptr<Primitive> aggregate,
        const std::vector<std::shared_ptr<Light>> &lights
    ) : lights(lights), aggregate(aggregate) {

        worldBound = aggregate->WorldBound();
    
        for (const auto &light : lights) {
            light->Preprocess(*this);
        }
    }

    const Bounds3f &WorldBound() const {
        return worldBound;
    }

    bool Intersect(const Ray &ray, SurfaceInteraction *si) const;

    // P is for "predicate". No intersection details are returned.
    bool IntersectP(const Ray &ray) const;

    // Tests the ray for intersection with the scene, including regions filled with
    // participating media. The beam transmittance at the point of intersection is
    // returned as well (the fraction of radiance that is transmitted between the
    // ray's origin and the point of intersection, which is subject to attenuation
    // due to absorption and out-scattering). 
    bool IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *si, Spectrum *transmittance) const;
};

#endif // CPBRT_CORE_SCENE_H