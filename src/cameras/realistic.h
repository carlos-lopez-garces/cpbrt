#ifndef CPBRT_CAMERAS_REALISTIC_H
#define CPBRT_CAMERAS_realistic_H

#include "core/cpbrt.h"
#include "core/camera.h"

class RealisticCamera : public Camera {
private:
    struct LensElementInterface {
        Float curvatureRadius;
        Float thickness;
        Float eta;
        Float apertureRadius;
    };

    const bool simpleWeighting;

    std::vector<LensElementInterface> elementInterfaces;

    Float LensRearZ() const;

    Float LensFrontZ() const;

    Float RearElementRadius() const;

    bool TraceLensesFromFilm(const Ray &ray, Ray *rayOut) const;

    bool TraceLensesFromScene(const Ray &rCamera, Ray *rayOut) const;

    void ComputeThickLensApproximation(Float pz[2], Float f[2]) const;

    Float FocusThickLens(Float focusDistance);

    Bounds2f BoundExitPupil(Float pFilmX0, Float pFilmX1) const;

    Point3f SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample, Float *sampleBoundsArea) const;

    Float FocusBinarySearch(Float focusDistance);

    static bool IntersectSphericalElement(
        Float radius,
        Float zCenter,
        const Ray &ray,
        Float *t,
        Normal3f *n
    );

    static void ComputeCardinalPoints(const Ray &rIn, const Ray &rOut, Float *p, Float *f);

public:
    RealisticCamera(
        const Transform &CameraToWorld,
        Float shutterOpen,
        Float shutterClose,
        Float apertureDiameter,
        Float focusDistance,
        bool simpleWeighting,
        std::vector<Float> &lensData,
        Film *film,
        const Medium *medium
    );

    Float GenerateRay(const CameraSample &sample, Ray *) const;
};

RealisticCamera *CreateRealisticCamera(
    const ParamSet &params, const Transform &cam2world, Film *film, const Medium *medium
);

#endif // CPBRT_CAMERAS_REALISTIC_H