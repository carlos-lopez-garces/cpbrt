#include "interaction.h"
#include "transform.h"
#include "primitive.h"
#include "shape.h"
#include "light.h"

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
        n = Faceforward(n, shading.n);
    } else {
        shading.n = Faceforward(shading.n, n);
    }

    shading.dpdu = dpdus;
    shading.dpdv = dpdvs;
    shading.dndu = dndus;
    shading.dndv = dndvs;
}

void SurfaceInteraction::ComputeDifferentials(const RayDifferential &ray) const {
    if (ray.hasDifferentials) {
        // Estimate screen space change in p and (u,v).

        // At the differential scale, the surface is assumed to be locally flat and can
        // be approximated by the plane ax + by + cz + d = 0 that is tangent to the surface
        // at p, the point of intersection.
        //
        // The coefficients of the plane implicit equation are known: a=nx, b=ny, c=nz, and
        // d=-dot(n,p), where n is the normal at p.
        //
        // The dpdx and dpdy differentials are computed as the difference vectors between p
        // and the intersection points of the ray differentials and the local approximation plane.
        // To compute these intersection points, the t parameters of the ray differentials are
        // computed as follows:
        Float d = -Dot(n, Vector3f(p.x, p.y, p.z));
        Float tx = (-Dot(n, Vector3f(ray.rxOrigin)) - d) / Dot(n, ray.rxDirection);
        Float ty = (-Dot(n, Vector3f(ray.ryOrigin)) - d) / Dot(n, ray.ryDirection);
        Point3f px = ray.rxOrigin + tx*ray.rxDirection;
        Point3f py = ray.ryOrigin + ty*ray.ryDirection;
        dpdx = px - p;
        dpdy = py - p;

        // Choose two dimensions to use for ray offset computation.
        // TODO: explain.
        int dim[2];
        if (std::abs(n.x) > std::abs(n.y) && std::abs(n.x) > std::abs(n.z)) {
            dim[0] = 1;
            dim[1] = 2;
        } else if (std::abs(n.y) > std::abs(n.z)) {
            dim[0] = 0;
            dim[1] = 2;
        } else {
            dim[0] = 0; 
            dim[1] = 1;
        }

        Float A[2][2] = {
            { dpdu[dim[0]], dpdv[dim[0]] },
            { dpdu[dim[1]], dpdv[dim[1]] }
        };

        Float bx[2] = {
            px[dim[0]] - p[dim[0]],
            px[dim[1]] - p[dim[1]]
        };

        Float by[2] = {
            py[dim[0]] - p[dim[0]],
            py[dim[1]] - p[dim[1]]
        };

        // The following 2 systems of equations are equivalent to this one:
        //
        // [ dp_x/du dp_x/dv ][ du ]   [ p'_x - p_x ]
        // [ dp_y/du dp_y/dv ][ dv ] = [ p'_y - p_y ]
        // [ dp_z/du dp_z/dv ]         [ p'_z - p_z ]
        //
        // which is the matrix equation form Ax=b of the vector equation 
        // p' - p = du dp/du + dv dp/dv that expresses either one of the intersection points
        // of the ray differentials and the local approximation plane as a linear combination
        // vector of the basis {dp/du, dpd/v} of UV parametric space.
        //
        // We are interested in finding (du/dx, dv/dx), the differential change in U and V
        // caused by a differential change in the px - p direction, where px is the ray x
        // differential intersection with the plane; and (du/dy, dv/dy) in the py - p direction,
        // where py is the intersection point of the ray y differential. 
        if (!SolveLinearSystem2x2(A, bx, &dudx, &dvdx)) {
            dudx = dvdx = 0;
        }
        if (!SolveLinearSystem2x2(A, by, &dudy, &dvdy)) {
            dudy = dvdy = 0;
        }
    } else {
        // 0-valued derivatives lead to unfiltered point sampling of textures.
        dudx = dudy = 0;
        dvdx = dvdy = 0;
        dpdx = dpdy = Vector3f(0, 0, 0);
    }
}

void SurfaceInteraction::ComputeScatteringFunctions(
    const RayDifferential &ray,
    MemoryArena &arena,
    bool allowMultipleLobes,
    TransportMode mode
) {
    // Compute differentials of P and UV mappings. They may be needed next by materials
    // that evaluate textures to obtain BSDF parameters.
    ComputeDifferentials(ray);

    // Delegate BSDF creation to primitive.
    primitive->ComputeScatteringFunctions(this, arena, mode, allowMultipleLobes);
}

Spectrum SurfaceInteraction::Le(const Vector3f &w) const {
    const AreaLight *areaLight = primitive->GetAreaLight();
    return areaLight ? areaLight->L(*this, w) : Spectrum(0.0f);
}