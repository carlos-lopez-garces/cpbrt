#include "orthographic.h"

Float OrthographicCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    // Compute the ray's origin by transforming the raster space sample's position
    // on the film to camera space.
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);

    // Rays generated by an orthographic camera always point down the z axis.
    *ray = Ray(pCamera, Vector3f(0, 0, 1));

    // TODO: Modify ray for depth of field.

    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;

    // Transform the ray to world space before returning it.
    *ray = CameraToWorld(*ray);

    // ?
    return 1;
}

Float OrthographicCamera::GenerateRayDifferential(
    const CameraSample &sample,
    RayDifferential *rd
) const {
    // Compute the main ray's origin by transforming the raster space sample's position
    // on the film to camera space.
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);

    // Rays generated by an orthographic camera always point down the z axis.
    *rd = RayDifferential(pCamera, Vector3f(0, 0, 1));

    // TODO: Modify ray for depth of field.

    if (lensRadius > 0) {
        // TODO: Compute OrthographicCamera ray differentials accounting for lens.
    } else {
        // The same camera space differentials apply to all rays and have been precomputed
        // by the constructor. 
        rd->rxOrigin = rd->o + dxCamera;
        rd->ryOrigin = rd->o + dyCamera;

        // Ray differentials are parallel to the main ray.
        rd->rxDirection = rd->d;
        rd->ryDirection = rd->d;
    }

    rd->time = Lerp(sample.time, shutterOpen, shutterClose);
    rd->hasDifferentials = true;
    rd->medium = medium;

    *rd = CameraToWorld(*rd);

    // ?
    return 1;
}