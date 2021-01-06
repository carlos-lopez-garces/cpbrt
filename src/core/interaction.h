#include "cpbrt.h"
#include "geometry.h"
#include "memory.h"
#include "transform.h"

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

    bool isSurfaceInteraction() const {
        return n != Normal3f();
    }
};

class SurfaceInteraction : public Interaction {
public:
    // UV coordinates from the parameterization of the surface.
    Point2f uv;

    // Partial derivatives of the point and the surface normal.
    Vector3f dpdu, dpdv;
    Normal3f dndu, dndv;

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

    // Initiates the creation of the BSDF at the surface-ray intersection point.
    void ComputeScatteringFunctions(
        const RayDifferential &ray,
        MemoryArena &arena,
        bool allowMultipleLobes,
        TransportMode mode
    );
};