#include "geometry.h"

template <typename T> inline bool Bounds3<T>::IntersectP(
    const Ray &ray,
    Float *hitt0,
    Float *hitt1
) const {
    // The parametric interval [t0, t1] shrinks with every slab intersection that is found.
    Float t0 = 0;
    Float t1 = ray.tMax;

    // Iteration 1 tests the intersection of the ray with the slab between the planes that
    // are perpendicular to the X axis. Iteration 2, Y axis. Iteration 3, Z axis.
    for (int i=0; i<3; ++i) {
        // Update interval for ith bounding box slab.

        // Avoid division later on. If component is 0, division by 0 of any nonzero floating-
        // point number results in a special +infinity or -infinity value in IEEE floating-
        // point arithmetic. The algorithm is supposed to be correct in the face of such case.
        Float reciprocalDirComponent = 1 / ray.d[i];
        
        // Compute values of the ray's parameter for which the ray satisfies the implicit
        // equation of the plane of. The ray intersects the slab at points (r.o + tNear*r.d)
        // and (r.o + tFar*r.d), if at all.
        Float tNear = (pMin[i] - ray.o[i]) *  reciprocalDirComponent;
        Float tFar = (pMax[i] - ray.o[i]) * reciprocalDirComponent;
        if (tNear > tFar) std::swap(tNear, tFar);

        // TODO: update tFar to ensure robust ray-bounds intersection (see 3.9.2).

        // Update (shrink) the parametric interval from slab intersection t values: 
        // [t0, t1] = [tNear, tFar].
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1) {
            // The parametric interval is empty; it shrunk until the endpoints crossed each other.
            return false;
        }
    }

    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}

// Optimized overload. 15% performance improvement approximately.
template <typename T> inline bool Bounds3<T>::IntersectP(
    const Ray &ray, 
    const Vector3f &reciprocalDir, 
    const int dirIsNeg[3]
) const {
    // TODO: implement.
    return false;
}
