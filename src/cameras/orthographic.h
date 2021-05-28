#ifndef CPBRT_CAMERAS_ORTHOGRAPHIC_H
#define CPBRT_CAMERAS_ORTHOGRAPHIC_H

#include "core/cpbrt.h"
#include "core/camera.h"
#include "core/film.h"

class OrthographicCamera : public ProjectiveCamera {
private:
    // x and y differentials used to shift the origin of ray differentials. Their
    // directions are parallel to the main ray's.
    Vector3f dxCamera;
    Vector3f dyCamera;

public:
    OrthographicCamera(
        const Transform &CameraToWorld,
        const Bounds2f &screenWindow,
        Float shutterOpen,
        Float shutterClose,
        Float lensRadius,
        Float focalDistance,
        Film *film,
        const Medium *medium
    ) : ProjectiveCamera(
            CameraToWorld,
            Orthographic(0, 1),
            screenWindow,
            shutterOpen,
            shutterClose,
            lensRadius,
            focalDistance,
            film,
            medium
    ) {
        // Compute differential changes in origin for orthographic camera rays. 
        // OrthographicCamera generates ray differentials that are parallel to the
        // main ray, their directions are the same; what it changes is their origin,
        // which it shifts always by the same x and y differentials: the camera space
        // distance that corresponds to 1 pixel in raster space.
        dxCamera = RasterToCamera(Vector3f(1, 0, 0));
        dyCamera = RasterToCamera(Vector3f(0, 1, 0));
    }

    // The orthographic transformation 1) translates camera space points along the z
    // axis so that the origin of camera space lies at the z coordinate of the near
    // plane (x and y coordinates remain unchanged); and 2) normalizes the z coordinate
    // of the points, between the near plane (x,y,0) and the far plane (x,y,1).
    Transform Orthographic(Float zNear, Float zFar) {
        return 
            Scale(1, 1, 1 / (zFar - zNear))
            * Translate(Vector3f(0, 0, -zNear));
    }

    Float GenerateRay(const CameraSample &sample, Ray *ray) const;

    Float GenerateRayDifferential(const CameraSample &sample, RayDifferential *rd) const;
};

#endif // CPBRT_CAMERAS_ORTHOGRAPHIC_H