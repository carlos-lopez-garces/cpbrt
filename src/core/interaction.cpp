#include "interaction.h"

SurfaceInteraction::SurfaceInteraction(
    const Point3f &p,
    const Vector3f &pError,
    const Point2f &uv,
    const Vector3f &wo,
    const Vector3f &dpdu,
    const Vector3f &dpdv,
    const Normal3f &dndu,
    const Normal3f &dndv,
    Float time,
    const Shape *shape
) : Interaction(
        p,
        // dpdu and dpdv lie on the tangent plane; their cross product is the normal, but it
        // might need to have its orientation and/or handedness flipped based on the Shape's
        // reverseOrientation attribute and on whether the Shape's Transform changed the 
        // handedness of the Shape's coordinate system.
        Normal3f(Normalize(Cross(dpdu, dpdv))),
        pError,
        wo,
        time,
        nullptr
), uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv), shape(shape) {
    shading.n = n;
    shading.dpdu = dpdu;
    shading.dpdv = dpdv;
    shading.dndu = dndu;
    shading.dndv = dndv;

    // Adjust normal based on orientation and handedness: the Shape may have been configured to
    // have its normal point inward or the Shape's Transform might change the coordinate system 
    // handedness.
    if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
        n *= -1;
        shading.n *= -1;
    }
}

// Set or update shading geometry. If the orientation of the new shading geometry is
// authoritative, it will override the orientation of the true geometry.
void SurfaceInteraction::SetShadingGeometry(
    const Vector3f &dpdus,
    const Vector3f &dpdvs,
    const Normal3f &dndus,
    const Normal3f &dndvs,
    bool orientationIsAuthoritative
) {
    shading.n = Normal3f(Normalize(Cross(dpdus, dpdvs)));
    if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
        shading.n *= -1;
    }

    if (orientationIsAuthoritative) {
        // Make the true geometric normal lie on the same hemisphere as the shading
        // normal.
        n = FaceForward(n, shading.n);
    } else {
        shading.n = FaceForward(shading.n, n)
    }

    shading.dpdu = dpdus;
    shading.dpdv = dpdvs;
    shading.dndu = dndus;
    shading.dndv = dndvs;
}

void SurfaceInteraction::ComputeScatteringFunctions(
    const RayDifferential &ray,
    MemoryArena &arena,
    bool allowMultipleLobes,
    TransportMode mode
) {
    // TODO: implement when implementing Texture.
    ComputeDifferentials(ray);

    // Delegate to Primitive.
    primitive->ComputeScatteringFunctions(this, arena, mode, allowMultipleLobes);
}