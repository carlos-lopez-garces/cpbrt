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

public:
    RealisticCamera(
        const AnimatedTransform &CameraToWorld,
        Float shutterOpen,
        Float shutterClose,
        Float apertureDiameter,
        Float focusDistance,
        bool simpleWeighting,
        std::vector<Float> &lensData,
        Film *film,
        const Medium *medium
    );
};

#endif // CPBRT_CAMERAS_REALISTIC_H