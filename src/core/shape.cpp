#include "shape.h"

Shape::Shape(
    const Transform *ObjectToWorld,
    const Transform *WorldToObject,
    bool reverseOrientation
) : ObjectToWorld(ObjectToWorld), 
    WorldToObject(WorldToObject), 
    reverseOrientation(reverseOrientation),
    transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
}

Float Shape::Pdf(const Interaction &it, const Vector3f &wi) const {
    // Does the sample ray from the given point (on some other surface) intersect
    // this Shape? This PDF is defined only over the set of directions along which
    // the Shape is visible from the input Interaction point.
    Ray ray = it.SpawnRay(wi);
    Float tHit;
    SurfaceInteraction shapeIntersection;
    if (!Intersect(ray, &tHit, &shapeIntersection, false)) {
        return 0;
    }

    // Note the 1/Area() factor, which corresponds to the PDF of a uniform distribution
    // of points over area (as computed by the other Shape::Pdf). Here we are
    // transforming that PDF over area into one over solid angle to sample a point that
    // is visible from the Interaction point.
    Float pdf = DistanceSquared(it.p, shapeIntersection.p) / (AbsDot(shapeIntersection.n, -wi) * Area());

    return pdf;
}