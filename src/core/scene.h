#ifndef CPBRT_CORE_SCENE_H
#define CPBRT_CORE_SCENE_H

#include <memory>

#include "geometry.h"
#include "light.h"
#include "primitive.h"

class Scene {
private:
    std::vector<std::shared_ptr<Light>> lights;

    // All the primitives in the scene.
    std::shared_ptr<Primitive> aggregate;

    Bounds3f worldBound;

public:
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
};

#endif // CPBRT_CORE_SCENE_H