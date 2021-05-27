#ifndef CPBRT_CORE_CAMERA_H
#define CPBRT_CORE_CAMERA_H

#include "cpbrt.h"
#include "film.h"
#include "geometry.h"
#include "transform.h"

struct CameraSample {
    // Point on the film.
    Point2f pFilm;

    // Point on the lens.
    Point2f pLens;

    // Point in time between the camera's [shutterOpen,shutterClose] interval at
    // which this sample is taken by a ray.
    Float time;
};

class Camera {
public:
    // TODO: The CameraToWorld transformation should really be an AnimatedTransform.
    Transform CameraToWorld;

    // Times at which the shutter opens and closes, respectively, defining an interval
    // of time during which the camera samples the scene.
    const Float shutterOpen;
    const Float shutterClose;

    // The captured image.
    Film *film;

    // The scattering medium in which the camera is immersed (typically air).
    const Medium *medium;

    Camera(
        const Transform &CameraToWorld,
        Float shutterOpen,
        Float shutterClose,
        Film *film,
        const Medium *medium
    ) : CameraToWorld(CameraToWorld),
        shutterOpen(shutterOpen),
        shutterClose(shutterClose),
        film(film),
        medium(medium) {
    }

    virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;

    Float GenerateRayDifferential(const CameraSample &sample, RayDifferential *rd) const;
};

class ProjectiveCamera : public Camera {
protected:
    // Implementations of ProjectiveCamera supply the CameraToScreen transformation.
    // The rest are computed by the ProjectiveCamera constructor.
    Transform CameraToScreen;
    Transform RasterToCamera;

    // Screen and raster space are isomorphic.
    Transform ScreenToRaster;
    Transform RasterToScreen;

    // Aperture. If 0, the camera behaves like in the pinhole model and every point is
    // imaged in focus regardless of its depth. If >0, only the points at depths in the
    // lens' depth of field are imaged in focus; the rest are blurred.
    Float lensRadius;

    // Camera space scene depth of plane of focus.
    Float focalDistance;
public:
    ProjectiveCamera(
        const Transform &CameraToWorld,
        const Transform &CameraToScreen,
        const Bounds2f  &screenWindow,
        Float shutterOpen,
        Float shutterClose,
        Float lensRadius,
        Float focalDistance,
        Film *film,
        const Medium *medium
    ) : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
        CameraToScreen(CameraToScreen),
        lensRadius(lensRadius),
        focalDistance(focalDistance) {
        // Compute projective camera transformations.

        // The ScreenToRaster transformation translates the screen space origin to 
        // the upper-left corner of the near plane raster (which corresponds to the
        // upper-left corner of the film), and along with it the rest of the vector
        // space.
        //
        // Then it normalizes the screen space coordinate system to obtain NDC space.
        //
        // And, finally, the transformation scales the NDC space up by the raster's
        // resolution to obtain raster space.
        Float screenWidth = screenWindow.pMax.x - screenWindow.pMin.x;
        // Note that the Y axis gets inverted.
        Float screenHeight = screenWindow.pMin.y - screenWindow.pMax.y;
        ScreenToRaster = 
            Scale(film->fullResolution.x, film->fullResolution.y, 1)
            * Scale(1 / screenWidth, 1 / screenHeight, 1)
            * Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));

        RasterToScreen = Inverse(ScreenToRaster);
        RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
    }
};

#endif // CPBRT_CORE_CAMERA_H