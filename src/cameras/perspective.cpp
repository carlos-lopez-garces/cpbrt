#include "perspective.h"
#include "core/paramset.h"
#include "core/sampler.h"
#include "core/sampling.h"
#include "core/light.h"

PerspectiveCamera::PerspectiveCamera(
    const Transform &CameraToWorld,
    const Bounds2f &screenWindow,
    Float shutterOpen,
    Float shutterClose,
    Float lensRadius,
    Float focalDistance,
    Float fov,
    Film *film,
    const Medium *medium
) : ProjectiveCamera(
    CameraToWorld,
    Perspective(fov, 1e-2f, 1000.f),
    screenWindow,
    shutterOpen,
    shutterClose,
    lensRadius,
    focalDistance,
    film,
    medium
) {
    // Compute differential changes in origin for perspective camera rays.
    dxCamera = RasterToCamera(Point3f(1, 0, 0)) - RasterToCamera(Point3f(0, 0, 0));
    dyCamera = RasterToCamera(Point3f(0, 1, 0)) - RasterToCamera(Point3f(0, 0, 0));

    // TODO: Compute image plane bounds at z=1 for PerspectiveCamera.
}

Float PerspectiveCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    // Compute the ray's direction by transforming the raster space sample's position
    // on the film to camera space. The vector passes through that point.
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);

    // Rays generated by a perspective camera have their origin at the origin of 
    // camera space.
    *ray = Ray(Point3f(0, 0, 0), Normalize(Vector3f(pCamera)));

    // Modify ray for effect of lens and depth of field. Note that this function operates
    // in camera space, not in screen or raster space, and its execution takes place
    // before the projection, so none of the effects of the perspective transformation
    // are relevant here. 
    if (lensRadius > 0) {
        // TODO: if the point being sampled is on the focal plane or in the DoF range,
        // does ConcentricSampleDisk with the input sample.pLens also return a jittered
        // output? Wouldn't that cause the point to be imaged out of focus?
        //
        // pLens is a randomized point on the lens that gives the ray its origin.
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // A ray uses the parametric definition of the line. We denote the ray's
        // parameter with t.
        //
        // ft is the ray's parameter that yields a point on the focal plane along the ray.
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Given a CameraSample, all the rays generated for it will point at the exact same
        // point of focus in the scene, but their origins will be randomized over the disk 
        // surface of the lens.
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;

    // Transform the ray to world space before returning it.
    *ray = CameraToWorld(*ray);

    // Cameras that model a system of lenses may assign different weights to different rays;
    // the weight is used to scale the ray's contribution to the image.
    // The PerspectiveCamera gives equal weight to all rays.
    return 1;
}

Float PerspectiveCamera::GenerateRayDifferential(
    const CameraSample &sample,
    RayDifferential *rd
) const {
    // Compute the main ray's direction by transforming the raster space sample's position
    // on the film to camera space. The vector passes through that point.
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);

    // Rays generated by a perspective camera have their origin at the origin of 
    // camera space.
    *rd = RayDifferential(Point3f(0, 0, 0), Normalize(Vector3f(pCamera)));

    // Modify main ray for effect of lens and depth of field. See comments in GenerateRay().
    if (lensRadius > 0) {
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        Float ft = focalDistance / rd->d.z;
        Point3f pFocus = (*rd)(ft);

        rd->o = Point3f(pLens.x, pLens.y, 0);
        rd->d = Normalize(pFocus - rd->o);
    }

    if (lensRadius > 0) {
        // Compute PerspectiveCamera ray differentials accounting for lens.
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Float ft = focalDistance / rd->d.z;

        Point3f pFocus = pCamera + dxCamera + (ft * Vector3f(0, 0, 1));
        rd->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        rd->rxDirection = Normalize(pFocus - rd->rxOrigin);

        pFocus = pCamera + dyCamera + (ft * Vector3f(0, 0, 1));
        rd->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        rd->ryDirection = Normalize(pFocus - rd->ryOrigin);
    } else {
        rd->rxOrigin = rd->o;
        rd->ryOrigin = rd->o;
        rd->rxDirection = Normalize(Vector3f(pCamera) + dxCamera);
        rd->ryDirection = Normalize(Vector3f(pCamera) + dyCamera);
    }

    rd->time = Lerp(sample.time, shutterOpen, shutterClose);
    rd->hasDifferentials = true;
    rd->medium = medium;

    // Transform the ray to world space before returning it.
    *rd = CameraToWorld(*rd);

    // Cameras that model a system of lenses may assign different weights to different rays;
    // the weight is used to scale the ray's contribution to the image.
    // The PerspectiveCamera gives equal weight to all rays.
    return 1;
}