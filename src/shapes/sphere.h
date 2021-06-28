#ifndef CPBRT_SHAPES_SPHERE_H
#define CPBRT_SHAPES_SPHERE_H

#include "core/cpbrt.h"
#include "core/geometry.h"
#include "core/shape.h"
#include "core/transform.h"

class Sphere : public Shape {
private:
    const Float radius;
    // Clipping plane z=zMin that clips the bottom of the sphere.
    const Float zMin;
    // Cipping plane z=Max that clips the top of the sphere.
    const Float zMax;
    // The Theta angle is measured from the xy-plane to the z-axis.
    // Only points (Theta, Phi) => (x, y, z) such that thetaMin < Theta < thetaMax are part of the
    // sphere.
    const Float thetaMin;
    const Float thetaMax;
    // Maximum polar angle. Only points (Theta, Phi) => (x, y, z) such that 0 < Phi < phiMax are part
    // of the sphere.
    const Float phiMax;

public:
    Sphere(
        const Transform *ObjectToWorld,
        const Transform *WorldToObject,
        bool reverseOrientation,
        Float radius,
        // The z=zMin and z=zMax planes that clip the sphere are used to determine the thetaMin
        // and thetaMax angles.
        Float zMin,
        Float zMax,
        Float phiMax
    ) : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        radius(radius),
        zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
        zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
        thetaMin(std::acos(Clamp(zMin / radius, -1, 1))),
        thetaMax(std::acos(Clamp(zMax / radius, -1, 1))),
        phiMax(Radians(Clamp(phiMax, 0, 360)))
    {}

    // Compute the object-space bounding box of the sphere.
    Bounds3f ObjectBound() const;

    bool Intersect(
        const Ray &r,
        Float *tHit,
        SurfaceInteraction *si,
        bool testAlphaTexture
    ) const;

    bool IntersectP(
        const Ray &r,
        bool testAlphaTexture
    ) const;

    Float Area() const;

    Interaction Sample(const Point2f &u, Float *pdf) const;

    Interaction Sample(const Interaction &ref, const Point2f &u, Float *pdf) const;

    Float Pdf(const Interaction &ref, const Vector3f &wi) const;
};

std::shared_ptr<Shape> CreateSphereShape(
    // Object to world.
    const Transform *o2w,
    // World to object.
    const Transform * w2o,
    bool reverseOrientation,
    const ParamSet &params
);

#endif // CPBRT_SHAPES_SPHERE_H