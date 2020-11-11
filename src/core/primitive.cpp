#include "primitive.h"

Bounds3f GeometricPrimitive::WorldBound() const {
    return shape->WorldBound();
}

bool GeometricPrimitive::Intersect(const Ray &ray, SurfaceInteraction &si) const {
    Float tHit;

    if (!shape->Intersect(ray, &tHit, si)) {
        return false;
    }

    // Subsequent Shapes that lie beyond tHit will be ignored in the remaining intersection tests.
    ray.tMax = tHit;
    si->primitive = this;
    // TODO: initialize SurfaceInteraction::mediumInterface after Shape intersection.
    return true;
}

// P is for "predicate". No intersection details are returned.
bool GeometricPrimitive::IntersectP(const Ray &ray) const {
    return shape->IntersectP(ray);
}

const AreaLight *GeometricPrimitive::GetAreaLight() const {
    return areaLight.get();
}

const Material *GeometricPrimitive::GetMaterial() const {
    return material.get();
}

void GeometricPrimitive::ComputeScatteringFunctions(
    SurfaceInteraction *si,
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    if (material) {
        material->ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    }
}