#ifndef CPBRT_CORE_SHAPE_H
#define CPBRT_CORE_SHAPE_H

#include "cpbrt.h"
#include "interaction.h"
#include "transform.h"
#include "geometry.h"

class Shape {
public:
    // WorldToObject is the inverse of ObjectToWorld. This pair of Transforms is an
    // isomorphism between object space and world space.
    // The memory of these Transforms is not managed by the Shape.
    const Transform *ObjectToWorld;
    const Transform *WorldToObject;

    // By default, normals point outwards. A Shape may be configured to have its normal point
    // inward.
    const bool reverseOrientation;
    const bool transformSwapsHandedness;

    Shape(
        const Transform *ObjectToWorld,
        const Transform *WorldToObject,
        bool reverseOrientation
    );

    // Bounding volume in object space.
    virtual Bounds3f ObjectBound() const = 0;

    // Bounding volume in world space.
    Bounds3f WorldBound() const;

    // Returns the details of the closest intersection between the ray and this Shape's geometry.
    // The input ray is always in world space and the returned intersection details should be in
    // world space as well.
    virtual bool Intersect(
        const Ray &ray,
        Float *tHit,
        SurfaceInteraction *isect,
        bool testAlphaTexture = true
    ) const = 0;

    // Tells whether an intersection occurs, without returning details. P is for "predicate".
    // This method should be an efficient way to tell if an intersection occurs. This default
    // implementation is not efficient because it calls Intersect, which computes the details
    // of the intersection. Classes implementing this interface should not do that, in general.
    virtual bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const {
        Float tHit = ray.tMax;
        SurfaceInteraction isect;
        return Intersect(ray, &tHit, &isect, testAlphaTexture);
    }

    virtual Float Area() const = 0;

    // Samples a point on the surface of the Shape using the input random sample value. The
    // sample is drawn according to a PDF defined with respect to the surface area of the Shape.
    virtual Interaction Sample(const Point2f &u) const = 0;

    // Evaluates a PDF defined with respect to area. A Shape has a uniform distribution over area
    // by default.
    virtual Float Pdf(const Interaction &) const {
        return 1 / Area();
    }

    // Samples the subset of points on the surface of the Shape that are visible from the input
    // Interaction point (which is a point on some other surface). The sample is drawn according
    // to a PDF defined with respect to solid angle (the sample is a point on the Shape, but the
    // corresponding probability density is with respect to the set of directions along which the
    // points are visible). 
    virtual Interaction Sample(const Interaction &it, const Point2f &u) const {
        return Sample(u);
    }

    // Evaluates a PDF defined with respect to solid angle: the set of directions along which the
    // shape is visible from the given Interaction point. 
    virtual Float Pdf(const Interaction &it, const Vector3f &wi) const;    
};

#endif // CPBRT_CORE_SHAPE_H