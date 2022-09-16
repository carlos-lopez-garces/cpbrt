#ifndef CPBRT_CAMERAS_REALISTIC_H
#define CPBRT_CAMERAS_realistic_H

#include "core/cpbrt.h"
#include "core/camera.h"

class RealisticCamera : public Camera {
private:

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