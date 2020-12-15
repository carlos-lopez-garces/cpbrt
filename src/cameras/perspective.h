#include <cmath>

#include "core/cpbrt.h"
#include "core/camera.h"

class PerspectiveCamera : public ProjectiveCamera {
private:
    Vector3f dxCamera;
    Vector3f dyCamera;

public:
    PerspectiveCamera(
        const Transform &CameraToWorld,
        const Bounds2f &screenWindow,
        Float shutterOpen,
        Float shutterClose,
        Float lensRadius,
        Float focalDistance,
        Float fov,
        Film *film,
        const Medium *medium
    );

    Transform Perspective(Float fov, Float zNear, Float zFar) {
        // Set up the canonical perpsective matrix.
        //
        // The x and y coordinates are mapped to x/z and y/z to create the 
        // foreshortening effect: the farther the point is, the larger its
        // camera space z coordinate is and the smaller the quotient will be,
        // drawing the point closer to the origin in screen space. The division
        // is achieved by the 4th column (0, 0, 1, 0), which causes the weight
        // entry of the homogeneous coordinates of the input vector to be z
        // after the multiplication.
        //
        // The z coordinate is normalized; the z coordinate of the near plane
        // being mapped to 0 and the far plane's to 1.
        Matrix4x4 perspective(
            1, 0, 0,                     0,
            0, 1, 0,                     0,
            0, 0, zFar / (zFar - zNear), -(zFar * zNear) / (zFar - zNear),
            0, 0, 1,                     0
        );

        // Scale canonical perspective view to specified field of view.
        // TODO: explain.
        Float reciprocalTanAng = 1 / std::tan(Radians(fov) / 2);

        return Scale(reciprocalTanAng, reciprocalTanAng, 1) * Transform(perspective);
    }

    Float GenerateRay(const CameraSample &sample, Ray *ray) const;

    Float GenerateRayDifferential(const CameraSample &sample, RayDifferential *rd) const;
};