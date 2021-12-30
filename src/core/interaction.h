#ifndef CPBRT_CORE_INTERACTION_H
#define CPBRT_CORE_INTERACTION_H

#include "cpbrt.h"
#include "geometry.h"
#include "transform.h"
#include "medium.h"
#include "material.h"
#include "error.h"

class Interaction {
public:
    Point3f p;

    Normal3f n;

    Float time;

    // An Interaction carries the floating-point error bound of the associated surface or medium.
    Vector3f pError;

    // Stands for 'omega sub o', the outgoing direction of light at the point p. The 0 vector
    // when 'omega sub o' doesn't apply for the Interaction's use case.
    Vector3f wo;

    // Scattering medium present at the point p.
    MediumInterface mediumInterface;

    Interaction() : time(0) {}

    Interaction(
        const Point3f &p,
        const Normal3f &n,
        const Vector3f &pError,
        const Vector3f &wo,
        Float time,
        const MediumInterface &mediumInterface
    ) : p(p), n(n), time(time), pError(pError), wo(wo), mediumInterface(mediumInterface) {}

    Interaction(
        const Point3f &p,
        Float time,
        const MediumInterface &mediumInterface
    ) : p(p), time(time), mediumInterface(mediumInterface) {}

    Interaction(
        const Point3f &p,
        const Vector3f &wo,
        Float time,
        const MediumInterface &mediumInterface
    ) : p(p), time(time), wo(wo), mediumInterface(mediumInterface) {}

    // Is the intersected primitive a surface?
    bool IsSurfaceInteraction() const {
        return n != Normal3f();
    }

    // Is the intersected primitive a region filled with participating media?
    bool IsMediumInteraction() const { 
        return !IsSurfaceInteraction();
    }

    Ray SpawnRay(const Vector3f &d) const {
        Point3f o = OffsetRayOrigin(p, pError, n, d);
        return Ray(o, d, Infinity, time, GetMedium(d));
    }

    Ray SpawnRayTo(const Point3f &p2) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
        Vector3f d = p2 - p;
        return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
    }

    Ray SpawnRayTo(const Interaction &it) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
        Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
        Vector3f d = target - origin;
        return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
    }

    // Returns the medium of the region where the ray will travel next after intersecting
    // this surface or region boundary.
    const Medium *GetMedium(const Vector3f &w) const {
        return Dot(w, n) > 0 ? mediumInterface.outside : mediumInterface.inside;
    }

    // Returns the medium where the intersection occurred, assuming that it occurred
    // inside a region filled with a participating medium, as opposed to at the
    // boundary of one or in a vacuum.
    const Medium *GetMedium() const {
        Assert(mediumInterface.inside == mediumInterface.outside);
        return mediumInterface.inside;
    }
};

class SurfaceInteraction : public Interaction {
public:
    // UV coordinates from the parameterization of the surface.
    Point2f uv;

    // Partial derivatives of the point and the surface normal.
    Vector3f dpdu, dpdv;
    Normal3f dndu, dndv;

    // Partial derivatives of the P: (x,y) -> (wpx, wpy, wpz) mapping from raster space
    // position to world space position with respect to raster space position, that is,
    // the change in the world space coordinate of p caused by a differential change in
    // the x and y directions in raster space.
    //
    // The change is computed using RayDifferentials, which are highly dependent on the
    // type of camera projection.
    mutable Vector3f dpdx, dpdy;

    // Partial derivatives of the U: (x,y) -> u and V: (x,y) -> v mappings from raster
    // space position to UV parametric space.
    mutable Float dudx = 0, dudy = 0;
    mutable Float dvdx = 0, dvdy = 0;

    // A Shape is assumed to have a parametric description: an isomorphism exists between
    // (a subspace of) R^2 and (a subspace of) R^3 that maps (u,v) coordinates onto
    // (x,y,z) points on the Shape.
    const Shape *shape = nullptr;

    const Primitive *primitive = nullptr;

    // Shading geometry. Perturbed normal and perturbed point and normal partial derivatives
    // are used by certain shading operations, such as bump mapping. The true normal is
    // Interaction::n and the true derivatives are elsewhere here in SurfaceInteraction.
    struct {
        Normal3f n;
        Vector3f dpdu;
        Vector3f dpdv;
        Normal3f dndu;
        Normal3f dndv;
    } shading;

    BSDF *bsdf = nullptr;
    BSSRDF *bssrdf = nullptr;

    SurfaceInteraction() {}

    SurfaceInteraction(
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
    );

    void SetShadingGeometry(
        const Vector3f &dpdus,
        const Vector3f &dpdvs,
        const Normal3f &dndus,
        const Normal3f &dndvs,
        bool orientationIsAuthoritative
    );

    // Computes the partial derivatives of the P: (x,y) -> (wpx, wpy, wpz) and UV mappings
    // U: (x,y) -> u and V: (x,y) -> v, where (x,y) is a point in raster space and
    // (wpx, wpy, wpz) is a point in world space.
    //
    // Ultimately what we want is to sample the texture function at the same rate as the
    // image function.
    void ComputeDifferentials(const RayDifferential &ray) const;

    // Initiates the creation of the BSDF at the surface-ray intersection point.
    void ComputeScatteringFunctions(
        const RayDifferential &ray,
        MemoryArena &arena,
        bool allowMultipleLobes = false,
        TransportMode mode = TransportMode::Radiance
    );

    // Computes the radiance emitted in the given direction by the intersected surface at
    // the intersection point. 0 if the surface is not an area light.
    Spectrum Le(const Vector3f &w) const;
};

class MediumInteraction : public Interaction {
public:
    // Describes the scattering (directional) distribution of the medium where
    // this intersection occurred.
    const PhaseFunction *phase;

    MediumInteraction(
        const Point3f &p,
        const Vector3f &wo,
        Float time,
        const Medium *medium,
        const PhaseFunction *phase
    ) : Interaction(p, wo, time, medium), phase(phase) {}

    bool IsValid() const {
        return phase != nullptr;
    }
};

#endif // CPBRT_CORE_INTERACTION_H