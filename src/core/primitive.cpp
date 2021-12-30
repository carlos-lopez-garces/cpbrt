#include "primitive.h"
#include "light.h"
#include "interaction.h"

Bounds3f GeometricPrimitive::WorldBound() const {
    return shape->WorldBound();
}

bool GeometricPrimitive::Intersect(const Ray &ray, SurfaceInteraction *si) const {
    Float tHit;

    if (!shape->Intersect(ray, &tHit, si)) {
        return false;
    }

    // Subsequent Shapes that lie beyond tHit will be ignored in the remaining intersection tests.
    ray.tMax = tHit;
    si->primitive = this;

    // Initialize SurfaceInteraction::mediumInterface after Shape intersection.
    if (mediumInterface.IsMediumTransition()) {
        // The primitive hit by the ray represents a region filled with a different scattering
        // medium and its surface is the boundary between 2 different scattering media. Record
        // the medium interface of the primitive in the hit.
        si->mediumInterface = mediumInterface;
    } else {
        // The ray didn't cross any boundary between scattering media. The ray origin
        // and the primitive are embedded in the same medium. 
        si->mediumInterface = MediumInterface(ray.medium);
    }

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

Bounds3f TransformedPrimitive::WorldBound() const {
    return primitive->WorldBound();
}

bool TransformedPrimitive::Intersect(const Ray &r, SurfaceInteraction *isect) const {
    // TODO: need to implement AnimatedTransform.
    return false;
}

bool TransformedPrimitive::IntersectP(const Ray &r) const {
    // TODO: need to implement AnimatedTransform.
    return false;
}

// Doesn't apply to an Aggregate. See SurfaceInteraction::primitive instead,
// where the Aggregate places the intersected Primitive.
const AreaLight *Aggregate::GetAreaLight() const {
    // TODO: log fatal error.
    return nullptr;
}

// Doesn't apply to an Aggregate. See SurfaceInteraction::primitive instead,
// where the Aggregate places the intersected Primitive.
const Material *Aggregate::GetMaterial() const {
    // TODO: log fatal error.
    return nullptr;
}

// Doesn't apply to an Aggregate. See SurfaceInteraction::primitive instead,
// where the Aggregate places the intersected Primitive.
void Aggregate::ComputeScatteringFunctions(
    SurfaceInteraction *si,
    MemoryArena &arena,
    TransportMode mode,
    bool allowMultipleLobes
) const {
    // TODO: log fatal error.
}