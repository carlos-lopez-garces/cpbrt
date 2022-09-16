#include "cameras/realistic.h"

RealisticCamera::RealisticCamera(
    const AnimatedTransform &CameraToWorld,
    Float shutterOpen,
    Float shutterClose,
    Float apertureDiameter,
    Float focusDistance,
    bool simpleWeighting,
    std::vector<Float> &lensData,
    Film *film,
    const Medium *medium
) : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
    simpleWeighting(simpleWeighting) 
{
    for (int i = 0; i < (int)lensData.size(); i += 4) {
        if (lensData[i] == 0) {
            if (apertureDiameter > lensData[i + 3]) {

            } else {
                lensData[i + 3] = apertureDiameter;
            }
        }

        elementInterfaces.push_back(
            LensElementInterface({lensData[i] * (Float).001, lensData[i + 1] * (Float).001, lensData[i + 2], lensData[i + 3] * Float(.001) / Float(2.)})    
        );
    }

    Float fb = FocusBinarySearch(focusDistance);
    elementInterfaces.back().thickness = FocusThickLens(focusDistance);

    int nSamples = 64;
    exitPupilBounds.resize(nSamples);
    ParallelFor([&](int i) {
        Float r0 = (Float)i / nSamples * film->diagonal / 2;
        Float r1 = (Float)(i + 1) / nSamples * film->diagonal / 2;
        exitPupilBounds[i] = BoundExitPupil(r0, r1);
    }, nSamples);
}

Float RealisticCamera::LensRearZ() const {
    return elementInterfaces.back().thickness;
}

Float RealisticCamera::LensFrontZ() const {
    Float zSum = 0;
    for (const LensElementInterface &element : elementInterfaces) {
        zSum += element.thickness;
    }
    return zSum;
}

Float RealisticCamera::RearElementRadius() const {
    return elementInterfaces.back().apertureRadius;
}

bool RealisticCamera::TraceLensesFromFilm(const Ray &rCamera, Ray *rayOut) const {
    Float elementZ = 0;

    static const Transform CameraToLens = Scale(1, 1, -1);

    Ray rLens = CameraToLens(rCamera);

    for (int i = elementInterfaces.size() - 1; i >= 0; --i) {
        const LensElementInterface &element = elementInterfaces[i];
        elementZ -= element.thickness;

        Float t;
        Normal3f n;

        bool isStop = (element.curvatureRadius == 0);
        if (isStop) {
            if (rLens.d.z >= 0.0) {
                return false;
            }
            t = (elementZ - rLens.o.z) / rLens.d.z;
        } else {
            Float radius = element.curvatureRadius;
            Float zCenter = elementZ + element.curvatureRadius;
            if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n)) {
                return false;
            }
        }

        Point3f pHit = rLens(t);
        Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
        if (r2 > element.apertureRadius * element.apertureRadius) {
            return false;
        }
        rLens.o = pHit;

        if (!isStop) {
            Vector3f w;
            Float etaI = element.eta;
            Float etaT = (i > 0 && elementInterfaces[i - 1].eta != 0) ? elementInterfaces[i - 1].eta : 1;
            if (!Refract(Normalize(-rLens.d), n, etaI / etaT, &w)) {
                return false;
            }
            rLens.d = w;
        }
    }

    if (rayOut != nullptr) {
        static const Transform LensToCamera = Scale(1, 1, -1);
        *rayOut = LensToCamera(rLens);
    }

    return true;
}