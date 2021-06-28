#include "core/efloat.h"
#include "core/paramset.h"
#include "sphere.h"

Bounds3f Sphere::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, zMin), Point3f(radius, radius, zMax));
}

bool Sphere::Intersect(
    const Ray &r,
    Float *tHit,
    SurfaceInteraction *si,
    bool testAlphaTexture
) const {
    Float phi;
    Point3f pHit;

    // Transform ray to object space. The transformation may introduce floating-point error
    // in the ray's origin and direction; the error is kept track of in oErr and dErr.
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Initialize EFloat ray coordinate values.
    // TODO: Implement EFloat.
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);

    // Compute quadratic sphere coefficients: at^2 + bt + c = 0. The equation results from
    // substituting the ray's parametric equation [x: ox + tdx, y: oy + tdy, z: oz + tdz] in
    // the sphere's implicit equation x^2 + y^2 + z^2 - r^2 = 0.
    EFloat a = dx*dx + dy*dy + dz*dz;
    EFloat b = 2 * (dx*ox + dy*oy + dz*oz);
    EFloat c = ox*ox + oy*oy + oz*oz - EFloat(radius)*EFloat(radius);

    // Solve quadratic equation for t values: up to 2 solutions may exist.
    EFloat t0, t1;
    // TODO: Implement Quadratic.
    // Quadratic may introduce floating-point error too.
    // Quadratic guarantees that t0 < t1 (if they exist).
    if (!Quadratic(a, b, c, &t0, &t1)) {
        // No solution, no intersection.
        return false;
    }

    // Check quadric shape t0 and t1 for nearest intersection. Assumes that t0 < t1.
    if (
        // If true, the closest intersection is beyond the ray's extent.
        t0.UpperBound() > ray.tMax
        // If true, the intersection or intersections are behind the ray's origin.
        || t1.LowerBound() <= 0
    ) {
        return false;
    }
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        // The closest intersection is behind the origin.
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) {
            // The farthest intersection is beyond the ray's extent.
            return false;
        }
    }

    // Compute sphere hit position and Phi. Ray's () operator computes the point along the ray
    // that corresponds to the value of the input t parameter. This point also lies on the surface
    // of the sphere.
    pHit = ray((Float) tShapeHit);
    // Refine sphere intersection point.
    // TODO: explain.
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) {
        // The hit is at one of the poles of the sphere. Displace the x coordinate a little so that
        // the arctan call doesn't divide by zero: Phi = arctan(y/x).
        pHit.x = 1e-5f * radius;
    }
    // Phi is the polar angle.
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) {
        // Map it to a positive angle.
        phi += 2 * Pi;
    }

    // Test sphere intersection against clipping parameters. The clipping parameters restrict the
    // z and Phi values of points of the surface of the sphere. If the t0 intersection point is not
    // part of the clipped surface, we can see if the t1 intersection is.
    if ((zMin > -radius && pHit.z < zMin)
        || (zMax < radius && pHit.z > zMax)
        || phi > phiMax) {

        if (tShapeHit == t1) {
            // The tested intersection was t1 already.
            return false;
        }
        if (t1.UpperBound() > ray.tMax) {
            return false;
        }
        tShapeHit = t1;

        // Compute sphere hit position and Phi, but now with t1. Note that this block of code is
        // exactly the same as above, which does it for t0.
        pHit = ray((Float) tShapeHit);
        // Refine sphere intersection point.
        // TODO: explain.
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) {
            // The hit is at one of the poles of the sphere. Displace the x coordinate a little so that
            // the arctan call doesn't divide by zero: Phi = arctan(y/x).
            pHit.x = 1e-5f * radius;
        }
        // Phi is the polar angle.
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) {
            // Map it to a positive angle.
            phi += 2 * Pi;
        }

        if ((zMin > -radius && pHit.z < zMin)
            || (zMax < radius && pHit.z > zMax)
            || phi > phiMax) {
            return false;
        }
    }

    // Find parametric representation of sphere hit. This is a uv parameterization of the surface
    // of the sphere that will be included in the returned SurfaceInteraction.

    // phi is guaranteed to be <= phiMax at this point. Normalize it / remap it to [0,1].
    Float u = phi / phiMax;

    // The theta angle is measured from the +z axis. Remap [-radius, radius] to [-1,1]
    // and normalize the z coordinate before computing the angle.
    Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Part of the parameterization are the partial derivatives of the point and the normal vector.
    
    // Compute partial derivatives of the position.
    // TODO: explain.
    Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    Float invZRadius = 1 / zRadius;
    Float cosPhi = pHit.x * invZRadius;
    Float sinPhi = pHit.y * invZRadius;
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv = (thetaMax - thetaMin) * Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

    // Compute partial derivatives of the normal. Used for antialiasing textures of object reflected
    // on spheres.
    // TODO: explain.
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv = (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) * Vector3f(pHit.x, pHit.y, pHit.z);
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv); 
    Float invEGF2 = 1 / (E * G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu + (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu + (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for sphere intersection.
    // TODO: explain.
    Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

    // Initialize SurfaceInteraction from parametric information.
    *si = (*ObjectToWorld)(SurfaceInteraction(
        pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, dndu, dndv, ray.time, this
    ));

    // Return the ray's parameter value of the intersection to the caller.
    *tHit = (Float) tShapeHit;

    return true;
}

bool Sphere::IntersectP(
    const Ray &r,
    bool testAlphaTexture
) const {
    Float phi;
    Point3f pHit;

    // Transform ray to object space. The transformation may introduce floating-point error
    // in the ray's origin and direction; the error is kept track of in oErr and dErr.
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Initialize EFloat ray coordinate values.
    // TODO: Implement EFloat.
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);

    // Compute quadratic sphere coefficients: at^2 + bt + c = 0. The equation results from
    // substituting the ray's parametric equation [x: ox + tdx, y: oy + tdy, z: oz + tdz] in
    // the sphere's implicit equation x^2 + y^2 + z^2 - r^2 = 0.
    EFloat a = dx*dx + dy*dy + dz*dz;
    EFloat b = 2 * (dx*ox + dy*oy + dz*oz);
    EFloat c = ox*ox + oy*oy + oz*oz - EFloat(radius)*EFloat(radius);

    // Solve quadratic equation for t values: up to 2 solutions may exist.
    EFloat t0, t1;
    // TODO: Implement Quadratic.
    // Quadratic may introduce floating-point error too.
    // Quadratic guarantees that t0 < t1 (if they exist).
    if (!Quadratic(a, b, c, &t0, &t1)) {
        // No solution, no intersection.
        return false;
    }

    // Check quadric shape t0 and t1 for nearest intersection. Assumes that t0 < t1.
    if (
        // If true, the closest intersection is beyond the ray's extent.
        t0.UpperBound() > ray.tMax
        // If true, the intersection or intersections are behind the ray's origin.
        || t1.LowerBound() <= 0
    ) {
        return false;
    }
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        // The closest intersection is behind the origin.
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) {
            // The farthest intersection is beyond the ray's extent.
            return false;
        }
    }

    // Compute sphere hit position and Phi. Ray's () operator computes the point along the ray
    // that corresponds to the value of the input t parameter. This point also lies on the surface
    // of the sphere.
    pHit = ray((Float) tShapeHit);
    // Refine sphere intersection point.
    // TODO: explain.
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) {
        // The hit is at one of the poles of the sphere. Displace the x coordinate a little so that
        // the arctan call doesn't divide by zero: Phi = arctan(y/x).
        pHit.x = 1e-5f * radius;
    }
    // Phi is the polar angle.
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) {
        // Map it to a positive angle.
        phi += 2 * Pi;
    }

    // Test sphere intersection against clipping parameters. The clipping parameters restrict the
    // z and Phi values of points of the surface of the sphere. If the t0 intersection point is not
    // part of the clipped surface, we can see if the t1 intersection is.
    if ((zMin > -radius && pHit.z < zMin)
        || (zMax < radius && pHit.z > zMax)
        || phi > phiMax) {

        if (tShapeHit == t1) {
            // The tested intersection was t1 already.
            return false;
        }
        if (t1.UpperBound() > ray.tMax) {
            return false;
        }
        tShapeHit = t1;

        // Compute sphere hit position and Phi, but now with t1. Note that this block of code is
        // exactly the same as above, which does it for t0.
        pHit = ray((Float) tShapeHit);
        // Refine sphere intersection point.
        // TODO: explain.
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) {
            // The hit is at one of the poles of the sphere. Displace the x coordinate a little so that
            // the arctan call doesn't divide by zero: Phi = arctan(y/x).
            pHit.x = 1e-5f * radius;
        }
        // Phi is the polar angle.
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) {
            // Map it to a positive angle.
            phi += 2 * Pi;
        }

        if ((zMin > -radius && pHit.z < zMin)
            || (zMax < radius && pHit.z > zMax)
            || phi > phiMax) {
            return false;
        }
    }

    // Find parametric representation of sphere hit. This is a uv parameterization of the surface
    // of the sphere that will be included in the returned SurfaceInteraction.

    // phi is guaranteed to be <= phiMax at this point. Normalize it / remap it to [0,1].
    Float u = phi / phiMax;

    // The theta angle is measured from the +z axis. Remap [-radius, radius] to [-1,1]
    // and normalize the z coordinate before computing the angle.
    Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Part of the parameterization are the partial derivatives of the point and the normal vector.
    
    // Compute partial derivatives of the position.
    // TODO: explain.
    Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    Float invZRadius = 1 / zRadius;
    Float cosPhi = pHit.x * invZRadius;
    Float sinPhi = pHit.y * invZRadius;
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv = (thetaMax - thetaMin) * Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

    // Compute partial derivatives of the normal. Used for antialiasing textures of object reflected
    // on spheres.
    // TODO: explain.
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv = (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) * Vector3f(pHit.x, pHit.y, pHit.z);
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv); 
    Float invEGF2 = 1 / (E * G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu + (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu + (f * F - g * E) * invEGF2 * dpdv);

    return true;
}

Float Sphere::Area() const {
    // This expression results from evaluating the definite integral that defines the area
    // of the surface of revolution that resuls from revolving the circular arc
    // f(z) = sqrt(r^2 - z^2) around the x axis.
    return phiMax * radius * (zMax - zMin);
}

std::shared_ptr<Shape> CreateSphereShape(
    // Object to world.
    const Transform *o2w,
    // World to object.
    const Transform * w2o,
    bool reverseOrientation,
    const ParamSet &params
) {
    Float radius = params.FindOneFloat("radius", 1.f);
    // If not set, don't clip the sphere.
    Float zmin = params.FindOneFloat("zmin", -radius);
    Float zmax = params.FindOneFloat("zmax", radius);
    Float phimax = params.FindOneFloat("phimax", 360.f);
    return std::make_shared<Sphere>(
        o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax
    );
}