#ifndef CPBRT_CORE_PRIMITIVE_H
#define CPBRT_CORE_PRIMITIVE_H

#include <memory>

#include "cpbrt.h"
#include "material.h"
#include "medium.h"
#include "transform.h"

class Primitive {
public:
    virtual Bounds3f WorldBound() const = 0;

    virtual bool Intersect(const Ray &ray, SurfaceInteraction &si) const = 0;
    // P is for "predicate". No intersection details are returned.
    virtual bool IntersectP(const Ray &ray) const = 0;

    virtual const AreaLight *GetAreaLight() const = 0;

    virtual const Material *GetMaterial() const = 0;

    // Includes in the SurfaceInteraction the BSDF and/or BSSRDF that describe scattering 
    // at the ray intersection point. The MemoryArena allocates memory for the BSDF and BSSRDF.
    virtual void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const = 0;
};

class GeometricPrimitive : public Primitive {
private:
    std::shared_ptr<Shape> shape;

    std::shared_ptr<Material> material;

    // Null when non-emissive.
    std::shared_ptr<AreaLight> areaLight;

    // Participating medium inside or outside the Primitive.
    MediumInterface mediumInterface;

public:
    GeometricPrimitive(
        const std::shared_ptr<Shape> &shape,
        const std::shared_ptr<Material> &material,
        const std::shared_ptr<AreaLight> &areaLight,
        const MediumInterface &mediumInterface)
    :   shape(shape), material(material), areaLight(areaLight), mediumInterface(mediumInterface) {}

    Bounds3f WorldBound() const;

    bool Intersect(const Ray &ray, SurfaceInteraction &si) const;
    // P is for "predicate". No intersection details are returned.
    bool IntersectP(const Ray &ray) const;

    // Creates the BSDF at the surface-ray intersection point.
    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

class TransformedPrimitive : public Primitive {
private:
    std::shared_ptr<Primitive> primitive;
    // Places the primitive in the scene, in world-space.
    // Should really be an AnimatedTransform, but I don't want to support animations.
    const Transform PrimitiveToWorld;

public:
    TransformedPrimitive(
        std::shared_ptr<Primitive> &primitive,
        const Transform &PrimitiveToWorld
    ) : primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {}

    bool Intersect(const Ray &ray, SurfaceInteraction *si) const;
};

class Aggregate : public Primitive {
public:
    // Doesn't apply. Shouldn't be called.
    const AreaLight *GetAreaLight() const;

    // Doesn't apply. Shouldn't be called.
    const Material *GetMaterial() const;

    // Doesn't apply. Shouldn't be called.
    void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const;
};

#endif // CPBRT_CORE_PRIMITIVE_H