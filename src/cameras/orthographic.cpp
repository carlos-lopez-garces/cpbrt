#include "orthographic.h"
#include "core/paramset.h"
#include "core/sampler.h"
#include "core/sampling.h"

Float OrthographicCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    // Compute the ray's origin by transforming the raster space sample's position
    // on the film to camera space.
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);

    // Rays generated by an orthographic camera always point down the z axis.
    *ray = Ray(pCamera, Vector3f(0, 0, 1));

    // Modify ray for effect of lens and depth of field. Note that this function operates
    // in camera space, not in screen or raster space, and its execution takes place before
    // the projection, so none of the effects of the orthographic transformation are
    // relevant here. 
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
        // Since z=1 for orthographic ray directions, ft=focalDistance.
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
    // The OrthographicCamera gives equal weight to all rays.
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

    // Modify main ray for effect of lens and depth of field. See comments in GenerateRay().
    if (lensRadius > 0) {
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        Float ft = focalDistance / rd->d.z;
        Point3f pFocus = (*rd)(ft);

        rd->o = Point3f(pLens.x, pLens.y, 0);
        rd->d = Normalize(pFocus - rd->o);
    }

    if (lensRadius > 0) {
        // Compute OrthographicCamera ray differentials accounting for lens.
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Float ft = focalDistance / rd->d.z;

        Point3f pFocus = pCamera + dxCamera + (ft * Vector3f(0, 0, 1));
        rd->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        rd->rxDirection = Normalize(pFocus - rd->rxOrigin);

        pFocus = pCamera + dyCamera + (ft * Vector3f(0, 0, 1));
        rd->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        rd->ryDirection = Normalize(pFocus - rd->ryOrigin);
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

    // Cameras that model a system of lenses may assign different weights to different rays;
    // the weight is used to scale the ray's contribution to the image.
    // The OrthographicCamera gives equal weight to all rays.
    return 1;
}

OrthographicCamera *CreateOrthographicCamera(const ParamSet &params,
                                             const Transform &cam2world,
                                             Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6f);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film->fullResolution.x) / Float(film->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    return new OrthographicCamera(cam2world, screen, shutteropen, shutterclose,
                                  lensradius, focaldistance, film, medium);
}