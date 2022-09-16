#include "cameras/realistic.h"

RealisticCamera::RealisticCamera(
    const Transform &CameraToWorld,
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

bool RealisticCamera::IntersectSphericalElement(
    Float radius,
    Float zCenter,
    const Ray &ray,
    Float *t,    
    Normal3f *n
) {
    Point3f o = ray.o - Vector3f(0, 0, zCenter);
    Float A = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
    Float B = 2 * (ray.d.x * o.x + ray.d.y * o.y + ray.d.z * o.z);
    Float C = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
    Float t0, t1;

    if (!Quadratic(A, B, C, &t0, &t1)) {
        return false;
    }

    bool useCloserT = (ray.d.z > 0) ^ (radius < 0);
    *t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
    if (*t < 0) {
        return false;
    }

    *n = Normal3f(Vector3f(o + *t * ray.d));
    *n = Faceforward(Normalize(*n), -ray.d);
    return true;
}

bool RealisticCamera::TraceLensesFromScene(const Ray &rCamera, Ray *rayOut) const {
    Float elementZ = -LensFrontZ();
    static const Transform CameraToLens = Scale(1, 1, -1);
    Ray rLens = CameraToLens(rCamera);

    for (size_t i = 0; i < elementInterfaces.size(); ++i) {
        const LensElementInterface &element = elementInterfaces[i];

        Float t;

        Normal3f n;
        bool isStop = (element.curvatureRadius == 0);
        if (isStop) {
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
            Vector3f wt;
            Float etaI = (i == 0 || elementInterfaces[i - 1].eta == 0) ? 1 : elementInterfaces[i - 1].eta;
            Float etaT = (elementInterfaces[i].eta != 0) ? elementInterfaces[i].eta : 1;
            if (!Refract(Normalize(-rLens.d), n, etaI / etaT, &wt)) {
                return false;
            }
            rLens.d = wt;
        }
        elementZ += element.thickness;
    }

    if (rayOut != nullptr) {
        static const Transform LensToCamera = Scale(1, 1, -1);
        *rayOut = LensToCamera(rLens);
    }

    return true;
}

void RealisticCamera::ComputeCardinalPoints(
    const Ray &rIn, const Ray &rayOut, Float *pz, Float *fz
) {
    Float tf = -rayOut.o.x / rayOut.d.x;
    *fz = -rayOut(tf).z;
    Float tp = (rIn.o.x - rayOut.o.x) / rayOut.d.x;
    *pz = -rayOut(tp).z;
}

void RealisticCamera::ComputeThickLensApproximation(Float pz[2], Float fz[2]) const {
    Float x = .001 * film->diagonal;

    Ray rScene(Point3f(x, 0, LensFrontZ() + 1), Vector3f(0, 0, -1));
    Ray rFilm;
    ComputeCardinalPoints(rScene, rFilm, &pz[0], &fz[0]);

    rFilm = Ray(Point3f(x, 0, LensRearZ() - 1), Vector3f(0, 0, 1));

    ComputeCardinalPoints(rFilm, rScene, &pz[1], &fz[1]);
}

Float RealisticCamera::FocusThickLens(Float focusDistance) {
    Float pz[2], fz[2];
    ComputeThickLensApproximation(pz, fz);
    Float f = fz[0] - pz[0];
    Float z = -focusDistance;
    Float c = (pz[1] - z - pz[0]) * (pz[1] - z - 4 * f - pz[0]);
    Float delta = 0.5f * (pz[1] - z + pz[0] - std::sqrt(c));
    return elementInterfaces.back().thickness + delta;
}

Bounds2f RealisticCamera::BoundExitPupil(Float pFilmX0, Float pFilmX1) const {
    Bounds2f pupilBounds;
    const int nSamples = 1024 * 1024;
    int nExitingRays = 0;

    Float rearRadius = RearElementRadius();
    Bounds2f projRearBounds(Point2f(-1.5f * rearRadius, -1.5f * rearRadius), Point2f(1.5f * rearRadius, 1.5f * rearRadius));

    for (int i = 0; i < nSamples; ++i) {
        Point3f pFilm(Lerp((i + 0.5f) / nSamples, pFilmX0, pFilmX1), 0, 0);
        Float u[2] = {RadicalInverse(0, i), RadicalInverse(1, i)};
        Point3f pRear(
            Lerp(u[0], projRearBounds.pMin.x, projRearBounds.pMax.x),
            Lerp(u[1], projRearBounds.pMin.y, projRearBounds.pMax.y),
            LensRearZ()
        );

        if (Inside(Point2f(pRear.x, pRear.y), pupilBounds) || TraceLensesFromFilm(Ray(pFilm, pRear - pFilm), nullptr)) {
            pupilBounds = Union(pupilBounds, Point2f(pRear.x, pRear.y));
            ++nExitingRays;
        }
    }

    if (nExitingRays == 0) {
        return projRearBounds;
    }

    pupilBounds = Expand(pupilBounds, 2 * projRearBounds.Diagonal().Length() / std::sqrt(nSamples));

    return pupilBounds;
}

Point3f RealisticCamera::SampleExitPupil(
    const Point2f &pFilm,
    const Point2f &lensSample,
    Float *sampleBoundsArea
) const {
    Float rFilm = std::sqrt(pFilm.x * pFilm.x + pFilm.y * pFilm.y);
    int rIndex = rFilm / (film->diagonal / 2) * exitPupilBounds.size();
    rIndex = std::min((int)exitPupilBounds.size() - 1, rIndex);
    Bounds2f pupilBounds = exitPupilBounds[rIndex];
    if (sampleBoundsArea) {
        *sampleBoundsArea = pupilBounds.Area();
    }

    Point2f pLens = pupilBounds.Lerp(lensSample);

    Float sinTheta = (rFilm != 0) ? pFilm.y / rFilm : 0;
    Float cosTheta = (rFilm != 0) ? pFilm.x / rFilm : 1;
    return Point3f(cosTheta * pLens.x - sinTheta * pLens.y, sinTheta * pLens.x + cosTheta * pLens.y, LensRearZ());
}

Float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    ++totalRays;

    Point2f s(sample.pFilm.x / film->fullResolution.x, sample.pFilm.y / film->fullResolution.y);
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, 0);

    Float exitPupilBoundsArea;
    Point3f pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), sample.pLens, &exitPupilBoundsArea);
    Ray rFilm(pFilm, pRear - pFilm, Infinity, Lerp(sample.time, shutterOpen, shutterClose));
    if (!TraceLensesFromFilm(rFilm, ray)) {
        ++vignettedRays;
        return 0;
    }

    *ray = CameraToWorld(*ray);
    ray->d = Normalize(ray->d);
    ray->medium = medium;

    Float cosTheta = Normalize(rFilm.d).z;
    Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
    if (simpleWeighting) {
        return cos4Theta * exitPupilBoundsArea / exitPupilBounds[0].Area();
    } else {
        return (shutterClose - shutterOpen) * (cos4Theta * exitPupilBoundsArea) / (LensRearZ() * LensRearZ());
    }
}

RealisticCamera *CreateRealisticCamera(
    const ParamSet &params,
    const Transform &cam2world,
    Film *film,
    const Medium *medium
) {
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        std::swap(shutterclose, shutteropen);
    }

    std::string lensFile = params.FindOneFilename("lensfile", "");
    Float apertureDiameter = params.FindOneFloat("aperturediameter", 1.0);
    Float focusDistance = params.FindOneFloat("focusdistance", 10.0);
    bool simpleWeighting = params.FindOneBool("simpleweighting", true);
    if (lensFile == "") {
        return nullptr;
    }

    std::vector<Float> lensData;
    if (!ReadFloatFile(lensFile.c_str(), &lensData)) {
        return nullptr;
    }

    if (lensData.size() % 4 != 0) {
        return nullptr;
    }

    return new RealisticCamera(
        cam2world, shutteropen, shutterclose, apertureDiameter, focusDistance, simpleWeighting, lensData, film, medium
    );
}

Float RealisticCamera::FocusBinarySearch(Float focusDistance) {
    Float filmDistanceLower;
    Film filmDistanceUpper;

    filmDistanceLower = filmDistanceUpper = FocusThickLens(focusDistance);
    while (FocusDistance(filmDistanceLower) > focusDistance) {
        filmDistanceLower *= 1.005f;
    }
    while (FocusDistance(filmDistanceUpper) < focusDistance) {
        filmDistanceUpper /= 1.005f;
    }

    for (int i = 0; i < 20; ++i) {
        Float fmid = 0.5f * (filmDistanceLower + filmDistanceUpper);
        Float midFocus = FocusDistance(fmid);
        if (midFocus < focusDistance) {
            filmDistanceLower = fmid;
        } else {
            filmDistanceUpper = fmid;
        }
    }

    return 0.5f * (filmDistanceLower + filmDistanceUpper);
}