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
    Bounds3f WorldBound() const {
        // Transforming a Bounds3f involves computing a new AABB that contains the original
        // AABB's corner points. What the world-space AABB actually bounds is then the
        // object-space AABB; as a result, the Shape may not be tightly bound by the 
        // world-space AABB.
        return (*ObjectToWorld)(ObjectBound());
    }
};