#include <array>

#include "cameras/realistic.h"
#include "core/efloat.h"
#include "core/floatfile.h"
#include "core/imageio.h"
#include "core/paramset.h"
#include "core/reflection.h"
#include "core/sampler.h"
#include "core/sampling.h"

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
    const Ray &rIn, const Ray &rOut, Float *pz, Float *fz
) {
    Float tf = -rOut.o.x / rOut.d.x;
    *fz = -rOut(tf).z;
    Float tp = (rIn.o.x - rOut.o.x) / rOut.d.x;
    *pz = -rOut(tp).z;
}

void RealisticCamera::ComputeThickLensApproximation(Float pz[2], Float fz[2]) const {
    Float x = .001 * film->diagonal;

    Ray rScene(Point3f(x, 0, LensFrontZ() + 1), Vector3f(0, 0, -1));
    Ray rFilm;
    TraceLensesFromScene(rScene, &rFilm);
    ComputeCardinalPoints(rScene, rFilm, &pz[0], &fz[0]);

    rFilm = Ray(Point3f(x, 0, LensRearZ() - 1), Vector3f(0, 0, 1));
    TraceLensesFromFilm(rFilm, &rScene);
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
    Point2f s(sample.pFilm.x / film->fullResolution.x, sample.pFilm.y / film->fullResolution.y);
    Point2f pFilm2 = film->GetPhysicalExtent().Lerp(s);
    Point3f pFilm(-pFilm2.x, pFilm2.y, 0);

    Float exitPupilBoundsArea;
    Point3f pRear = SampleExitPupil(Point2f(pFilm.x, pFilm.y), sample.pLens, &exitPupilBoundsArea);
    Ray rFilm(pFilm, pRear - pFilm, Infinity, Lerp(sample.time, shutterOpen, shutterClose));
    if (!TraceLensesFromFilm(rFilm, ray)) {
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
    Float filmDistanceUpper;

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

Float RealisticCamera::FocusDistance(Float filmDistance) {
    Bounds2f bounds = BoundExitPupil(0, .001 * film->diagonal);

    const std::array<Float, 3> scaleFactors = {0.1f, 0.01f, 0.001f};
    Float lu = 0.0f;

    Ray ray;

    bool foundFocusRay = false;
    for (Float scale : scaleFactors) {
        lu = scale * bounds.pMax[0];
        if (TraceLensesFromFilm(Ray(Point3f(0, 0, LensRearZ() - filmDistance), Vector3f(lu, 0, filmDistance)), &ray)) {
            foundFocusRay = true;
            break;
        }
    }

    if (!foundFocusRay) {
        return Infinity;
    }

    Float tFocus = -ray.o.x / ray.d.x;
    Float zFocus = ray(tFocus).z;
    if (zFocus < 0) {
        zFocus = Infinity;
    }
    return zFocus;
}

inline uint32_t ReverseBits32(uint32_t n) {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    return n;
}

inline uint64_t ReverseBits64(uint64_t n) {
    uint64_t n0 = ReverseBits32((uint32_t)n);
    uint64_t n1 = ReverseBits32((uint32_t)(n >> 32));
    return (n0 << 32) | n1;
}

template <int base> static Float RadicalInverseSpecialized(uint64_t a) {
    const Float invBase = (Float)1 / (Float)base;
    uint64_t reversedDigits = 0;
    Float invBaseN = 1;
    while (a) {
        uint64_t next = a / base;
        uint64_t digit = a - next * base;
        reversedDigits = reversedDigits * base + digit;
        invBaseN *= invBase;
        a = next;
    }
    return std::min(reversedDigits * invBaseN, OneMinusEpsilon);
}

Float RadicalInverse(int baseIndex, uint64_t a) {
    switch (baseIndex) {
    case 0:
// #ifndef CPBRT_HAVE_HEX_FP_CONSTANTS
//         return ReverseBits64(a) * 5.4210108624275222e-20;
// #else
        return ReverseBits64(a) * 0x1p-64;
// #endif
    case 1:
        return RadicalInverseSpecialized<3>(a);
    case 2:
        return RadicalInverseSpecialized<5>(a);
    case 3:
        return RadicalInverseSpecialized<7>(a);
    // Remainder of cases for _RadicalInverse()_
    case 4:
        return RadicalInverseSpecialized<11>(a);
    case 5:
        return RadicalInverseSpecialized<13>(a);
    case 6:
        return RadicalInverseSpecialized<17>(a);
    case 7:
        return RadicalInverseSpecialized<19>(a);
    case 8:
        return RadicalInverseSpecialized<23>(a);
    case 9:
        return RadicalInverseSpecialized<29>(a);
    case 10:
        return RadicalInverseSpecialized<31>(a);
    case 11:
        return RadicalInverseSpecialized<37>(a);
    case 12:
        return RadicalInverseSpecialized<41>(a);
    case 13:
        return RadicalInverseSpecialized<43>(a);
    case 14:
        return RadicalInverseSpecialized<47>(a);
    case 15:
        return RadicalInverseSpecialized<53>(a);
    case 16:
        return RadicalInverseSpecialized<59>(a);
    case 17:
        return RadicalInverseSpecialized<61>(a);
    case 18:
        return RadicalInverseSpecialized<67>(a);
    case 19:
        return RadicalInverseSpecialized<71>(a);
    case 20:
        return RadicalInverseSpecialized<73>(a);
    case 21:
        return RadicalInverseSpecialized<79>(a);
    case 22:
        return RadicalInverseSpecialized<83>(a);
    case 23:
        return RadicalInverseSpecialized<89>(a);
    case 24:
        return RadicalInverseSpecialized<97>(a);
    case 25:
        return RadicalInverseSpecialized<101>(a);
    case 26:
        return RadicalInverseSpecialized<103>(a);
    case 27:
        return RadicalInverseSpecialized<107>(a);
    case 28:
        return RadicalInverseSpecialized<109>(a);
    case 29:
        return RadicalInverseSpecialized<113>(a);
    case 30:
        return RadicalInverseSpecialized<127>(a);
    case 31:
        return RadicalInverseSpecialized<131>(a);
    case 32:
        return RadicalInverseSpecialized<137>(a);
    case 33:
        return RadicalInverseSpecialized<139>(a);
    case 34:
        return RadicalInverseSpecialized<149>(a);
    case 35:
        return RadicalInverseSpecialized<151>(a);
    case 36:
        return RadicalInverseSpecialized<157>(a);
    case 37:
        return RadicalInverseSpecialized<163>(a);
    case 38:
        return RadicalInverseSpecialized<167>(a);
    case 39:
        return RadicalInverseSpecialized<173>(a);
    case 40:
        return RadicalInverseSpecialized<179>(a);
    case 41:
        return RadicalInverseSpecialized<181>(a);
    case 42:
        return RadicalInverseSpecialized<191>(a);
    case 43:
        return RadicalInverseSpecialized<193>(a);
    case 44:
        return RadicalInverseSpecialized<197>(a);
    case 45:
        return RadicalInverseSpecialized<199>(a);
    case 46:
        return RadicalInverseSpecialized<211>(a);
    case 47:
        return RadicalInverseSpecialized<223>(a);
    case 48:
        return RadicalInverseSpecialized<227>(a);
    case 49:
        return RadicalInverseSpecialized<229>(a);
    case 50:
        return RadicalInverseSpecialized<233>(a);
    case 51:
        return RadicalInverseSpecialized<239>(a);
    case 52:
        return RadicalInverseSpecialized<241>(a);
    case 53:
        return RadicalInverseSpecialized<251>(a);
    case 54:
        return RadicalInverseSpecialized<257>(a);
    case 55:
        return RadicalInverseSpecialized<263>(a);
    case 56:
        return RadicalInverseSpecialized<269>(a);
    case 57:
        return RadicalInverseSpecialized<271>(a);
    case 58:
        return RadicalInverseSpecialized<277>(a);
    case 59:
        return RadicalInverseSpecialized<281>(a);
    case 60:
        return RadicalInverseSpecialized<283>(a);
    case 61:
        return RadicalInverseSpecialized<293>(a);
    case 62:
        return RadicalInverseSpecialized<307>(a);
    case 63:
        return RadicalInverseSpecialized<311>(a);
    case 64:
        return RadicalInverseSpecialized<313>(a);
    case 65:
        return RadicalInverseSpecialized<317>(a);
    case 66:
        return RadicalInverseSpecialized<331>(a);
    case 67:
        return RadicalInverseSpecialized<337>(a);
    case 68:
        return RadicalInverseSpecialized<347>(a);
    case 69:
        return RadicalInverseSpecialized<349>(a);
    case 70:
        return RadicalInverseSpecialized<353>(a);
    case 71:
        return RadicalInverseSpecialized<359>(a);
    case 72:
        return RadicalInverseSpecialized<367>(a);
    case 73:
        return RadicalInverseSpecialized<373>(a);
    case 74:
        return RadicalInverseSpecialized<379>(a);
    case 75:
        return RadicalInverseSpecialized<383>(a);
    case 76:
        return RadicalInverseSpecialized<389>(a);
    case 77:
        return RadicalInverseSpecialized<397>(a);
    case 78:
        return RadicalInverseSpecialized<401>(a);
    case 79:
        return RadicalInverseSpecialized<409>(a);
    case 80:
        return RadicalInverseSpecialized<419>(a);
    case 81:
        return RadicalInverseSpecialized<421>(a);
    case 82:
        return RadicalInverseSpecialized<431>(a);
    case 83:
        return RadicalInverseSpecialized<433>(a);
    case 84:
        return RadicalInverseSpecialized<439>(a);
    case 85:
        return RadicalInverseSpecialized<443>(a);
    case 86:
        return RadicalInverseSpecialized<449>(a);
    case 87:
        return RadicalInverseSpecialized<457>(a);
    case 88:
        return RadicalInverseSpecialized<461>(a);
    case 89:
        return RadicalInverseSpecialized<463>(a);
    case 90:
        return RadicalInverseSpecialized<467>(a);
    case 91:
        return RadicalInverseSpecialized<479>(a);
    case 92:
        return RadicalInverseSpecialized<487>(a);
    case 93:
        return RadicalInverseSpecialized<491>(a);
    case 94:
        return RadicalInverseSpecialized<499>(a);
    case 95:
        return RadicalInverseSpecialized<503>(a);
    case 96:
        return RadicalInverseSpecialized<509>(a);
    case 97:
        return RadicalInverseSpecialized<521>(a);
    case 98:
        return RadicalInverseSpecialized<523>(a);
    case 99:
        return RadicalInverseSpecialized<541>(a);
    case 100:
        return RadicalInverseSpecialized<547>(a);
    case 101:
        return RadicalInverseSpecialized<557>(a);
    case 102:
        return RadicalInverseSpecialized<563>(a);
    case 103:
        return RadicalInverseSpecialized<569>(a);
    case 104:
        return RadicalInverseSpecialized<571>(a);
    case 105:
        return RadicalInverseSpecialized<577>(a);
    case 106:
        return RadicalInverseSpecialized<587>(a);
    case 107:
        return RadicalInverseSpecialized<593>(a);
    case 108:
        return RadicalInverseSpecialized<599>(a);
    case 109:
        return RadicalInverseSpecialized<601>(a);
    case 110:
        return RadicalInverseSpecialized<607>(a);
    case 111:
        return RadicalInverseSpecialized<613>(a);
    case 112:
        return RadicalInverseSpecialized<617>(a);
    case 113:
        return RadicalInverseSpecialized<619>(a);
    case 114:
        return RadicalInverseSpecialized<631>(a);
    case 115:
        return RadicalInverseSpecialized<641>(a);
    case 116:
        return RadicalInverseSpecialized<643>(a);
    case 117:
        return RadicalInverseSpecialized<647>(a);
    case 118:
        return RadicalInverseSpecialized<653>(a);
    case 119:
        return RadicalInverseSpecialized<659>(a);
    case 120:
        return RadicalInverseSpecialized<661>(a);
    case 121:
        return RadicalInverseSpecialized<673>(a);
    case 122:
        return RadicalInverseSpecialized<677>(a);
    case 123:
        return RadicalInverseSpecialized<683>(a);
    case 124:
        return RadicalInverseSpecialized<691>(a);
    case 125:
        return RadicalInverseSpecialized<701>(a);
    case 126:
        return RadicalInverseSpecialized<709>(a);
    case 127:
        return RadicalInverseSpecialized<719>(a);
    case 128:
        return RadicalInverseSpecialized<727>(a);
    case 129:
        return RadicalInverseSpecialized<733>(a);
    case 130:
        return RadicalInverseSpecialized<739>(a);
    case 131:
        return RadicalInverseSpecialized<743>(a);
    case 132:
        return RadicalInverseSpecialized<751>(a);
    case 133:
        return RadicalInverseSpecialized<757>(a);
    case 134:
        return RadicalInverseSpecialized<761>(a);
    case 135:
        return RadicalInverseSpecialized<769>(a);
    case 136:
        return RadicalInverseSpecialized<773>(a);
    case 137:
        return RadicalInverseSpecialized<787>(a);
    case 138:
        return RadicalInverseSpecialized<797>(a);
    case 139:
        return RadicalInverseSpecialized<809>(a);
    case 140:
        return RadicalInverseSpecialized<811>(a);
    case 141:
        return RadicalInverseSpecialized<821>(a);
    case 142:
        return RadicalInverseSpecialized<823>(a);
    case 143:
        return RadicalInverseSpecialized<827>(a);
    case 144:
        return RadicalInverseSpecialized<829>(a);
    case 145:
        return RadicalInverseSpecialized<839>(a);
    case 146:
        return RadicalInverseSpecialized<853>(a);
    case 147:
        return RadicalInverseSpecialized<857>(a);
    case 148:
        return RadicalInverseSpecialized<859>(a);
    case 149:
        return RadicalInverseSpecialized<863>(a);
    case 150:
        return RadicalInverseSpecialized<877>(a);
    case 151:
        return RadicalInverseSpecialized<881>(a);
    case 152:
        return RadicalInverseSpecialized<883>(a);
    case 153:
        return RadicalInverseSpecialized<887>(a);
    case 154:
        return RadicalInverseSpecialized<907>(a);
    case 155:
        return RadicalInverseSpecialized<911>(a);
    case 156:
        return RadicalInverseSpecialized<919>(a);
    case 157:
        return RadicalInverseSpecialized<929>(a);
    case 158:
        return RadicalInverseSpecialized<937>(a);
    case 159:
        return RadicalInverseSpecialized<941>(a);
    case 160:
        return RadicalInverseSpecialized<947>(a);
    case 161:
        return RadicalInverseSpecialized<953>(a);
    case 162:
        return RadicalInverseSpecialized<967>(a);
    case 163:
        return RadicalInverseSpecialized<971>(a);
    case 164:
        return RadicalInverseSpecialized<977>(a);
    case 165:
        return RadicalInverseSpecialized<983>(a);
    case 166:
        return RadicalInverseSpecialized<991>(a);
    case 167:
        return RadicalInverseSpecialized<997>(a);
    case 168:
        return RadicalInverseSpecialized<1009>(a);
    case 169:
        return RadicalInverseSpecialized<1013>(a);
    case 170:
        return RadicalInverseSpecialized<1019>(a);
    case 171:
        return RadicalInverseSpecialized<1021>(a);
    case 172:
        return RadicalInverseSpecialized<1031>(a);
    case 173:
        return RadicalInverseSpecialized<1033>(a);
    case 174:
        return RadicalInverseSpecialized<1039>(a);
    case 175:
        return RadicalInverseSpecialized<1049>(a);
    case 176:
        return RadicalInverseSpecialized<1051>(a);
    case 177:
        return RadicalInverseSpecialized<1061>(a);
    case 178:
        return RadicalInverseSpecialized<1063>(a);
    case 179:
        return RadicalInverseSpecialized<1069>(a);
    case 180:
        return RadicalInverseSpecialized<1087>(a);
    case 181:
        return RadicalInverseSpecialized<1091>(a);
    case 182:
        return RadicalInverseSpecialized<1093>(a);
    case 183:
        return RadicalInverseSpecialized<1097>(a);
    case 184:
        return RadicalInverseSpecialized<1103>(a);
    case 185:
        return RadicalInverseSpecialized<1109>(a);
    case 186:
        return RadicalInverseSpecialized<1117>(a);
    case 187:
        return RadicalInverseSpecialized<1123>(a);
    case 188:
        return RadicalInverseSpecialized<1129>(a);
    case 189:
        return RadicalInverseSpecialized<1151>(a);
    case 190:
        return RadicalInverseSpecialized<1153>(a);
    case 191:
        return RadicalInverseSpecialized<1163>(a);
    case 192:
        return RadicalInverseSpecialized<1171>(a);
    case 193:
        return RadicalInverseSpecialized<1181>(a);
    case 194:
        return RadicalInverseSpecialized<1187>(a);
    case 195:
        return RadicalInverseSpecialized<1193>(a);
    case 196:
        return RadicalInverseSpecialized<1201>(a);
    case 197:
        return RadicalInverseSpecialized<1213>(a);
    case 198:
        return RadicalInverseSpecialized<1217>(a);
    case 199:
        return RadicalInverseSpecialized<1223>(a);
    case 200:
        return RadicalInverseSpecialized<1229>(a);
    case 201:
        return RadicalInverseSpecialized<1231>(a);
    case 202:
        return RadicalInverseSpecialized<1237>(a);
    case 203:
        return RadicalInverseSpecialized<1249>(a);
    case 204:
        return RadicalInverseSpecialized<1259>(a);
    case 205:
        return RadicalInverseSpecialized<1277>(a);
    case 206:
        return RadicalInverseSpecialized<1279>(a);
    case 207:
        return RadicalInverseSpecialized<1283>(a);
    case 208:
        return RadicalInverseSpecialized<1289>(a);
    case 209:
        return RadicalInverseSpecialized<1291>(a);
    case 210:
        return RadicalInverseSpecialized<1297>(a);
    case 211:
        return RadicalInverseSpecialized<1301>(a);
    case 212:
        return RadicalInverseSpecialized<1303>(a);
    case 213:
        return RadicalInverseSpecialized<1307>(a);
    case 214:
        return RadicalInverseSpecialized<1319>(a);
    case 215:
        return RadicalInverseSpecialized<1321>(a);
    case 216:
        return RadicalInverseSpecialized<1327>(a);
    case 217:
        return RadicalInverseSpecialized<1361>(a);
    case 218:
        return RadicalInverseSpecialized<1367>(a);
    case 219:
        return RadicalInverseSpecialized<1373>(a);
    case 220:
        return RadicalInverseSpecialized<1381>(a);
    case 221:
        return RadicalInverseSpecialized<1399>(a);
    case 222:
        return RadicalInverseSpecialized<1409>(a);
    case 223:
        return RadicalInverseSpecialized<1423>(a);
    case 224:
        return RadicalInverseSpecialized<1427>(a);
    case 225:
        return RadicalInverseSpecialized<1429>(a);
    case 226:
        return RadicalInverseSpecialized<1433>(a);
    case 227:
        return RadicalInverseSpecialized<1439>(a);
    case 228:
        return RadicalInverseSpecialized<1447>(a);
    case 229:
        return RadicalInverseSpecialized<1451>(a);
    case 230:
        return RadicalInverseSpecialized<1453>(a);
    case 231:
        return RadicalInverseSpecialized<1459>(a);
    case 232:
        return RadicalInverseSpecialized<1471>(a);
    case 233:
        return RadicalInverseSpecialized<1481>(a);
    case 234:
        return RadicalInverseSpecialized<1483>(a);
    case 235:
        return RadicalInverseSpecialized<1487>(a);
    case 236:
        return RadicalInverseSpecialized<1489>(a);
    case 237:
        return RadicalInverseSpecialized<1493>(a);
    case 238:
        return RadicalInverseSpecialized<1499>(a);
    case 239:
        return RadicalInverseSpecialized<1511>(a);
    case 240:
        return RadicalInverseSpecialized<1523>(a);
    case 241:
        return RadicalInverseSpecialized<1531>(a);
    case 242:
        return RadicalInverseSpecialized<1543>(a);
    case 243:
        return RadicalInverseSpecialized<1549>(a);
    case 244:
        return RadicalInverseSpecialized<1553>(a);
    case 245:
        return RadicalInverseSpecialized<1559>(a);
    case 246:
        return RadicalInverseSpecialized<1567>(a);
    case 247:
        return RadicalInverseSpecialized<1571>(a);
    case 248:
        return RadicalInverseSpecialized<1579>(a);
    case 249:
        return RadicalInverseSpecialized<1583>(a);
    case 250:
        return RadicalInverseSpecialized<1597>(a);
    case 251:
        return RadicalInverseSpecialized<1601>(a);
    case 252:
        return RadicalInverseSpecialized<1607>(a);
    case 253:
        return RadicalInverseSpecialized<1609>(a);
    case 254:
        return RadicalInverseSpecialized<1613>(a);
    case 255:
        return RadicalInverseSpecialized<1619>(a);
    case 256:
        return RadicalInverseSpecialized<1621>(a);
    case 257:
        return RadicalInverseSpecialized<1627>(a);
    case 258:
        return RadicalInverseSpecialized<1637>(a);
    case 259:
        return RadicalInverseSpecialized<1657>(a);
    case 260:
        return RadicalInverseSpecialized<1663>(a);
    case 261:
        return RadicalInverseSpecialized<1667>(a);
    case 262:
        return RadicalInverseSpecialized<1669>(a);
    case 263:
        return RadicalInverseSpecialized<1693>(a);
    case 264:
        return RadicalInverseSpecialized<1697>(a);
    case 265:
        return RadicalInverseSpecialized<1699>(a);
    case 266:
        return RadicalInverseSpecialized<1709>(a);
    case 267:
        return RadicalInverseSpecialized<1721>(a);
    case 268:
        return RadicalInverseSpecialized<1723>(a);
    case 269:
        return RadicalInverseSpecialized<1733>(a);
    case 270:
        return RadicalInverseSpecialized<1741>(a);
    case 271:
        return RadicalInverseSpecialized<1747>(a);
    case 272:
        return RadicalInverseSpecialized<1753>(a);
    case 273:
        return RadicalInverseSpecialized<1759>(a);
    case 274:
        return RadicalInverseSpecialized<1777>(a);
    case 275:
        return RadicalInverseSpecialized<1783>(a);
    case 276:
        return RadicalInverseSpecialized<1787>(a);
    case 277:
        return RadicalInverseSpecialized<1789>(a);
    case 278:
        return RadicalInverseSpecialized<1801>(a);
    case 279:
        return RadicalInverseSpecialized<1811>(a);
    case 280:
        return RadicalInverseSpecialized<1823>(a);
    case 281:
        return RadicalInverseSpecialized<1831>(a);
    case 282:
        return RadicalInverseSpecialized<1847>(a);
    case 283:
        return RadicalInverseSpecialized<1861>(a);
    case 284:
        return RadicalInverseSpecialized<1867>(a);
    case 285:
        return RadicalInverseSpecialized<1871>(a);
    case 286:
        return RadicalInverseSpecialized<1873>(a);
    case 287:
        return RadicalInverseSpecialized<1877>(a);
    case 288:
        return RadicalInverseSpecialized<1879>(a);
    case 289:
        return RadicalInverseSpecialized<1889>(a);
    case 290:
        return RadicalInverseSpecialized<1901>(a);
    case 291:
        return RadicalInverseSpecialized<1907>(a);
    case 292:
        return RadicalInverseSpecialized<1913>(a);
    case 293:
        return RadicalInverseSpecialized<1931>(a);
    case 294:
        return RadicalInverseSpecialized<1933>(a);
    case 295:
        return RadicalInverseSpecialized<1949>(a);
    case 296:
        return RadicalInverseSpecialized<1951>(a);
    case 297:
        return RadicalInverseSpecialized<1973>(a);
    case 298:
        return RadicalInverseSpecialized<1979>(a);
    case 299:
        return RadicalInverseSpecialized<1987>(a);
    case 300:
        return RadicalInverseSpecialized<1993>(a);
    case 301:
        return RadicalInverseSpecialized<1997>(a);
    case 302:
        return RadicalInverseSpecialized<1999>(a);
    case 303:
        return RadicalInverseSpecialized<2003>(a);
    case 304:
        return RadicalInverseSpecialized<2011>(a);
    case 305:
        return RadicalInverseSpecialized<2017>(a);
    case 306:
        return RadicalInverseSpecialized<2027>(a);
    case 307:
        return RadicalInverseSpecialized<2029>(a);
    case 308:
        return RadicalInverseSpecialized<2039>(a);
    case 309:
        return RadicalInverseSpecialized<2053>(a);
    case 310:
        return RadicalInverseSpecialized<2063>(a);
    case 311:
        return RadicalInverseSpecialized<2069>(a);
    case 312:
        return RadicalInverseSpecialized<2081>(a);
    case 313:
        return RadicalInverseSpecialized<2083>(a);
    case 314:
        return RadicalInverseSpecialized<2087>(a);
    case 315:
        return RadicalInverseSpecialized<2089>(a);
    case 316:
        return RadicalInverseSpecialized<2099>(a);
    case 317:
        return RadicalInverseSpecialized<2111>(a);
    case 318:
        return RadicalInverseSpecialized<2113>(a);
    case 319:
        return RadicalInverseSpecialized<2129>(a);
    case 320:
        return RadicalInverseSpecialized<2131>(a);
    case 321:
        return RadicalInverseSpecialized<2137>(a);
    case 322:
        return RadicalInverseSpecialized<2141>(a);
    case 323:
        return RadicalInverseSpecialized<2143>(a);
    case 324:
        return RadicalInverseSpecialized<2153>(a);
    case 325:
        return RadicalInverseSpecialized<2161>(a);
    case 326:
        return RadicalInverseSpecialized<2179>(a);
    case 327:
        return RadicalInverseSpecialized<2203>(a);
    case 328:
        return RadicalInverseSpecialized<2207>(a);
    case 329:
        return RadicalInverseSpecialized<2213>(a);
    case 330:
        return RadicalInverseSpecialized<2221>(a);
    case 331:
        return RadicalInverseSpecialized<2237>(a);
    case 332:
        return RadicalInverseSpecialized<2239>(a);
    case 333:
        return RadicalInverseSpecialized<2243>(a);
    case 334:
        return RadicalInverseSpecialized<2251>(a);
    case 335:
        return RadicalInverseSpecialized<2267>(a);
    case 336:
        return RadicalInverseSpecialized<2269>(a);
    case 337:
        return RadicalInverseSpecialized<2273>(a);
    case 338:
        return RadicalInverseSpecialized<2281>(a);
    case 339:
        return RadicalInverseSpecialized<2287>(a);
    case 340:
        return RadicalInverseSpecialized<2293>(a);
    case 341:
        return RadicalInverseSpecialized<2297>(a);
    case 342:
        return RadicalInverseSpecialized<2309>(a);
    case 343:
        return RadicalInverseSpecialized<2311>(a);
    case 344:
        return RadicalInverseSpecialized<2333>(a);
    case 345:
        return RadicalInverseSpecialized<2339>(a);
    case 346:
        return RadicalInverseSpecialized<2341>(a);
    case 347:
        return RadicalInverseSpecialized<2347>(a);
    case 348:
        return RadicalInverseSpecialized<2351>(a);
    case 349:
        return RadicalInverseSpecialized<2357>(a);
    case 350:
        return RadicalInverseSpecialized<2371>(a);
    case 351:
        return RadicalInverseSpecialized<2377>(a);
    case 352:
        return RadicalInverseSpecialized<2381>(a);
    case 353:
        return RadicalInverseSpecialized<2383>(a);
    case 354:
        return RadicalInverseSpecialized<2389>(a);
    case 355:
        return RadicalInverseSpecialized<2393>(a);
    case 356:
        return RadicalInverseSpecialized<2399>(a);
    case 357:
        return RadicalInverseSpecialized<2411>(a);
    case 358:
        return RadicalInverseSpecialized<2417>(a);
    case 359:
        return RadicalInverseSpecialized<2423>(a);
    case 360:
        return RadicalInverseSpecialized<2437>(a);
    case 361:
        return RadicalInverseSpecialized<2441>(a);
    case 362:
        return RadicalInverseSpecialized<2447>(a);
    case 363:
        return RadicalInverseSpecialized<2459>(a);
    case 364:
        return RadicalInverseSpecialized<2467>(a);
    case 365:
        return RadicalInverseSpecialized<2473>(a);
    case 366:
        return RadicalInverseSpecialized<2477>(a);
    case 367:
        return RadicalInverseSpecialized<2503>(a);
    case 368:
        return RadicalInverseSpecialized<2521>(a);
    case 369:
        return RadicalInverseSpecialized<2531>(a);
    case 370:
        return RadicalInverseSpecialized<2539>(a);
    case 371:
        return RadicalInverseSpecialized<2543>(a);
    case 372:
        return RadicalInverseSpecialized<2549>(a);
    case 373:
        return RadicalInverseSpecialized<2551>(a);
    case 374:
        return RadicalInverseSpecialized<2557>(a);
    case 375:
        return RadicalInverseSpecialized<2579>(a);
    case 376:
        return RadicalInverseSpecialized<2591>(a);
    case 377:
        return RadicalInverseSpecialized<2593>(a);
    case 378:
        return RadicalInverseSpecialized<2609>(a);
    case 379:
        return RadicalInverseSpecialized<2617>(a);
    case 380:
        return RadicalInverseSpecialized<2621>(a);
    case 381:
        return RadicalInverseSpecialized<2633>(a);
    case 382:
        return RadicalInverseSpecialized<2647>(a);
    case 383:
        return RadicalInverseSpecialized<2657>(a);
    case 384:
        return RadicalInverseSpecialized<2659>(a);
    case 385:
        return RadicalInverseSpecialized<2663>(a);
    case 386:
        return RadicalInverseSpecialized<2671>(a);
    case 387:
        return RadicalInverseSpecialized<2677>(a);
    case 388:
        return RadicalInverseSpecialized<2683>(a);
    case 389:
        return RadicalInverseSpecialized<2687>(a);
    case 390:
        return RadicalInverseSpecialized<2689>(a);
    case 391:
        return RadicalInverseSpecialized<2693>(a);
    case 392:
        return RadicalInverseSpecialized<2699>(a);
    case 393:
        return RadicalInverseSpecialized<2707>(a);
    case 394:
        return RadicalInverseSpecialized<2711>(a);
    case 395:
        return RadicalInverseSpecialized<2713>(a);
    case 396:
        return RadicalInverseSpecialized<2719>(a);
    case 397:
        return RadicalInverseSpecialized<2729>(a);
    case 398:
        return RadicalInverseSpecialized<2731>(a);
    case 399:
        return RadicalInverseSpecialized<2741>(a);
    case 400:
        return RadicalInverseSpecialized<2749>(a);
    case 401:
        return RadicalInverseSpecialized<2753>(a);
    case 402:
        return RadicalInverseSpecialized<2767>(a);
    case 403:
        return RadicalInverseSpecialized<2777>(a);
    case 404:
        return RadicalInverseSpecialized<2789>(a);
    case 405:
        return RadicalInverseSpecialized<2791>(a);
    case 406:
        return RadicalInverseSpecialized<2797>(a);
    case 407:
        return RadicalInverseSpecialized<2801>(a);
    case 408:
        return RadicalInverseSpecialized<2803>(a);
    case 409:
        return RadicalInverseSpecialized<2819>(a);
    case 410:
        return RadicalInverseSpecialized<2833>(a);
    case 411:
        return RadicalInverseSpecialized<2837>(a);
    case 412:
        return RadicalInverseSpecialized<2843>(a);
    case 413:
        return RadicalInverseSpecialized<2851>(a);
    case 414:
        return RadicalInverseSpecialized<2857>(a);
    case 415:
        return RadicalInverseSpecialized<2861>(a);
    case 416:
        return RadicalInverseSpecialized<2879>(a);
    case 417:
        return RadicalInverseSpecialized<2887>(a);
    case 418:
        return RadicalInverseSpecialized<2897>(a);
    case 419:
        return RadicalInverseSpecialized<2903>(a);
    case 420:
        return RadicalInverseSpecialized<2909>(a);
    case 421:
        return RadicalInverseSpecialized<2917>(a);
    case 422:
        return RadicalInverseSpecialized<2927>(a);
    case 423:
        return RadicalInverseSpecialized<2939>(a);
    case 424:
        return RadicalInverseSpecialized<2953>(a);
    case 425:
        return RadicalInverseSpecialized<2957>(a);
    case 426:
        return RadicalInverseSpecialized<2963>(a);
    case 427:
        return RadicalInverseSpecialized<2969>(a);
    case 428:
        return RadicalInverseSpecialized<2971>(a);
    case 429:
        return RadicalInverseSpecialized<2999>(a);
    case 430:
        return RadicalInverseSpecialized<3001>(a);
    case 431:
        return RadicalInverseSpecialized<3011>(a);
    case 432:
        return RadicalInverseSpecialized<3019>(a);
    case 433:
        return RadicalInverseSpecialized<3023>(a);
    case 434:
        return RadicalInverseSpecialized<3037>(a);
    case 435:
        return RadicalInverseSpecialized<3041>(a);
    case 436:
        return RadicalInverseSpecialized<3049>(a);
    case 437:
        return RadicalInverseSpecialized<3061>(a);
    case 438:
        return RadicalInverseSpecialized<3067>(a);
    case 439:
        return RadicalInverseSpecialized<3079>(a);
    case 440:
        return RadicalInverseSpecialized<3083>(a);
    case 441:
        return RadicalInverseSpecialized<3089>(a);
    case 442:
        return RadicalInverseSpecialized<3109>(a);
    case 443:
        return RadicalInverseSpecialized<3119>(a);
    case 444:
        return RadicalInverseSpecialized<3121>(a);
    case 445:
        return RadicalInverseSpecialized<3137>(a);
    case 446:
        return RadicalInverseSpecialized<3163>(a);
    case 447:
        return RadicalInverseSpecialized<3167>(a);
    case 448:
        return RadicalInverseSpecialized<3169>(a);
    case 449:
        return RadicalInverseSpecialized<3181>(a);
    case 450:
        return RadicalInverseSpecialized<3187>(a);
    case 451:
        return RadicalInverseSpecialized<3191>(a);
    case 452:
        return RadicalInverseSpecialized<3203>(a);
    case 453:
        return RadicalInverseSpecialized<3209>(a);
    case 454:
        return RadicalInverseSpecialized<3217>(a);
    case 455:
        return RadicalInverseSpecialized<3221>(a);
    case 456:
        return RadicalInverseSpecialized<3229>(a);
    case 457:
        return RadicalInverseSpecialized<3251>(a);
    case 458:
        return RadicalInverseSpecialized<3253>(a);
    case 459:
        return RadicalInverseSpecialized<3257>(a);
    case 460:
        return RadicalInverseSpecialized<3259>(a);
    case 461:
        return RadicalInverseSpecialized<3271>(a);
    case 462:
        return RadicalInverseSpecialized<3299>(a);
    case 463:
        return RadicalInverseSpecialized<3301>(a);
    case 464:
        return RadicalInverseSpecialized<3307>(a);
    case 465:
        return RadicalInverseSpecialized<3313>(a);
    case 466:
        return RadicalInverseSpecialized<3319>(a);
    case 467:
        return RadicalInverseSpecialized<3323>(a);
    case 468:
        return RadicalInverseSpecialized<3329>(a);
    case 469:
        return RadicalInverseSpecialized<3331>(a);
    case 470:
        return RadicalInverseSpecialized<3343>(a);
    case 471:
        return RadicalInverseSpecialized<3347>(a);
    case 472:
        return RadicalInverseSpecialized<3359>(a);
    case 473:
        return RadicalInverseSpecialized<3361>(a);
    case 474:
        return RadicalInverseSpecialized<3371>(a);
    case 475:
        return RadicalInverseSpecialized<3373>(a);
    case 476:
        return RadicalInverseSpecialized<3389>(a);
    case 477:
        return RadicalInverseSpecialized<3391>(a);
    case 478:
        return RadicalInverseSpecialized<3407>(a);
    case 479:
        return RadicalInverseSpecialized<3413>(a);
    case 480:
        return RadicalInverseSpecialized<3433>(a);
    case 481:
        return RadicalInverseSpecialized<3449>(a);
    case 482:
        return RadicalInverseSpecialized<3457>(a);
    case 483:
        return RadicalInverseSpecialized<3461>(a);
    case 484:
        return RadicalInverseSpecialized<3463>(a);
    case 485:
        return RadicalInverseSpecialized<3467>(a);
    case 486:
        return RadicalInverseSpecialized<3469>(a);
    case 487:
        return RadicalInverseSpecialized<3491>(a);
    case 488:
        return RadicalInverseSpecialized<3499>(a);
    case 489:
        return RadicalInverseSpecialized<3511>(a);
    case 490:
        return RadicalInverseSpecialized<3517>(a);
    case 491:
        return RadicalInverseSpecialized<3527>(a);
    case 492:
        return RadicalInverseSpecialized<3529>(a);
    case 493:
        return RadicalInverseSpecialized<3533>(a);
    case 494:
        return RadicalInverseSpecialized<3539>(a);
    case 495:
        return RadicalInverseSpecialized<3541>(a);
    case 496:
        return RadicalInverseSpecialized<3547>(a);
    case 497:
        return RadicalInverseSpecialized<3557>(a);
    case 498:
        return RadicalInverseSpecialized<3559>(a);
    case 499:
        return RadicalInverseSpecialized<3571>(a);
    case 500:
        return RadicalInverseSpecialized<3581>(a);
    case 501:
        return RadicalInverseSpecialized<3583>(a);
    case 502:
        return RadicalInverseSpecialized<3593>(a);
    case 503:
        return RadicalInverseSpecialized<3607>(a);
    case 504:
        return RadicalInverseSpecialized<3613>(a);
    case 505:
        return RadicalInverseSpecialized<3617>(a);
    case 506:
        return RadicalInverseSpecialized<3623>(a);
    case 507:
        return RadicalInverseSpecialized<3631>(a);
    case 508:
        return RadicalInverseSpecialized<3637>(a);
    case 509:
        return RadicalInverseSpecialized<3643>(a);
    case 510:
        return RadicalInverseSpecialized<3659>(a);
    case 511:
        return RadicalInverseSpecialized<3671>(a);
    case 512:
        return RadicalInverseSpecialized<3673>(a);
    case 513:
        return RadicalInverseSpecialized<3677>(a);
    case 514:
        return RadicalInverseSpecialized<3691>(a);
    case 515:
        return RadicalInverseSpecialized<3697>(a);
    case 516:
        return RadicalInverseSpecialized<3701>(a);
    case 517:
        return RadicalInverseSpecialized<3709>(a);
    case 518:
        return RadicalInverseSpecialized<3719>(a);
    case 519:
        return RadicalInverseSpecialized<3727>(a);
    case 520:
        return RadicalInverseSpecialized<3733>(a);
    case 521:
        return RadicalInverseSpecialized<3739>(a);
    case 522:
        return RadicalInverseSpecialized<3761>(a);
    case 523:
        return RadicalInverseSpecialized<3767>(a);
    case 524:
        return RadicalInverseSpecialized<3769>(a);
    case 525:
        return RadicalInverseSpecialized<3779>(a);
    case 526:
        return RadicalInverseSpecialized<3793>(a);
    case 527:
        return RadicalInverseSpecialized<3797>(a);
    case 528:
        return RadicalInverseSpecialized<3803>(a);
    case 529:
        return RadicalInverseSpecialized<3821>(a);
    case 530:
        return RadicalInverseSpecialized<3823>(a);
    case 531:
        return RadicalInverseSpecialized<3833>(a);
    case 532:
        return RadicalInverseSpecialized<3847>(a);
    case 533:
        return RadicalInverseSpecialized<3851>(a);
    case 534:
        return RadicalInverseSpecialized<3853>(a);
    case 535:
        return RadicalInverseSpecialized<3863>(a);
    case 536:
        return RadicalInverseSpecialized<3877>(a);
    case 537:
        return RadicalInverseSpecialized<3881>(a);
    case 538:
        return RadicalInverseSpecialized<3889>(a);
    case 539:
        return RadicalInverseSpecialized<3907>(a);
    case 540:
        return RadicalInverseSpecialized<3911>(a);
    case 541:
        return RadicalInverseSpecialized<3917>(a);
    case 542:
        return RadicalInverseSpecialized<3919>(a);
    case 543:
        return RadicalInverseSpecialized<3923>(a);
    case 544:
        return RadicalInverseSpecialized<3929>(a);
    case 545:
        return RadicalInverseSpecialized<3931>(a);
    case 546:
        return RadicalInverseSpecialized<3943>(a);
    case 547:
        return RadicalInverseSpecialized<3947>(a);
    case 548:
        return RadicalInverseSpecialized<3967>(a);
    case 549:
        return RadicalInverseSpecialized<3989>(a);
    case 550:
        return RadicalInverseSpecialized<4001>(a);
    case 551:
        return RadicalInverseSpecialized<4003>(a);
    case 552:
        return RadicalInverseSpecialized<4007>(a);
    case 553:
        return RadicalInverseSpecialized<4013>(a);
    case 554:
        return RadicalInverseSpecialized<4019>(a);
    case 555:
        return RadicalInverseSpecialized<4021>(a);
    case 556:
        return RadicalInverseSpecialized<4027>(a);
    case 557:
        return RadicalInverseSpecialized<4049>(a);
    case 558:
        return RadicalInverseSpecialized<4051>(a);
    case 559:
        return RadicalInverseSpecialized<4057>(a);
    case 560:
        return RadicalInverseSpecialized<4073>(a);
    case 561:
        return RadicalInverseSpecialized<4079>(a);
    case 562:
        return RadicalInverseSpecialized<4091>(a);
    case 563:
        return RadicalInverseSpecialized<4093>(a);
    case 564:
        return RadicalInverseSpecialized<4099>(a);
    case 565:
        return RadicalInverseSpecialized<4111>(a);
    case 566:
        return RadicalInverseSpecialized<4127>(a);
    case 567:
        return RadicalInverseSpecialized<4129>(a);
    case 568:
        return RadicalInverseSpecialized<4133>(a);
    case 569:
        return RadicalInverseSpecialized<4139>(a);
    case 570:
        return RadicalInverseSpecialized<4153>(a);
    case 571:
        return RadicalInverseSpecialized<4157>(a);
    case 572:
        return RadicalInverseSpecialized<4159>(a);
    case 573:
        return RadicalInverseSpecialized<4177>(a);
    case 574:
        return RadicalInverseSpecialized<4201>(a);
    case 575:
        return RadicalInverseSpecialized<4211>(a);
    case 576:
        return RadicalInverseSpecialized<4217>(a);
    case 577:
        return RadicalInverseSpecialized<4219>(a);
    case 578:
        return RadicalInverseSpecialized<4229>(a);
    case 579:
        return RadicalInverseSpecialized<4231>(a);
    case 580:
        return RadicalInverseSpecialized<4241>(a);
    case 581:
        return RadicalInverseSpecialized<4243>(a);
    case 582:
        return RadicalInverseSpecialized<4253>(a);
    case 583:
        return RadicalInverseSpecialized<4259>(a);
    case 584:
        return RadicalInverseSpecialized<4261>(a);
    case 585:
        return RadicalInverseSpecialized<4271>(a);
    case 586:
        return RadicalInverseSpecialized<4273>(a);
    case 587:
        return RadicalInverseSpecialized<4283>(a);
    case 588:
        return RadicalInverseSpecialized<4289>(a);
    case 589:
        return RadicalInverseSpecialized<4297>(a);
    case 590:
        return RadicalInverseSpecialized<4327>(a);
    case 591:
        return RadicalInverseSpecialized<4337>(a);
    case 592:
        return RadicalInverseSpecialized<4339>(a);
    case 593:
        return RadicalInverseSpecialized<4349>(a);
    case 594:
        return RadicalInverseSpecialized<4357>(a);
    case 595:
        return RadicalInverseSpecialized<4363>(a);
    case 596:
        return RadicalInverseSpecialized<4373>(a);
    case 597:
        return RadicalInverseSpecialized<4391>(a);
    case 598:
        return RadicalInverseSpecialized<4397>(a);
    case 599:
        return RadicalInverseSpecialized<4409>(a);
    case 600:
        return RadicalInverseSpecialized<4421>(a);
    case 601:
        return RadicalInverseSpecialized<4423>(a);
    case 602:
        return RadicalInverseSpecialized<4441>(a);
    case 603:
        return RadicalInverseSpecialized<4447>(a);
    case 604:
        return RadicalInverseSpecialized<4451>(a);
    case 605:
        return RadicalInverseSpecialized<4457>(a);
    case 606:
        return RadicalInverseSpecialized<4463>(a);
    case 607:
        return RadicalInverseSpecialized<4481>(a);
    case 608:
        return RadicalInverseSpecialized<4483>(a);
    case 609:
        return RadicalInverseSpecialized<4493>(a);
    case 610:
        return RadicalInverseSpecialized<4507>(a);
    case 611:
        return RadicalInverseSpecialized<4513>(a);
    case 612:
        return RadicalInverseSpecialized<4517>(a);
    case 613:
        return RadicalInverseSpecialized<4519>(a);
    case 614:
        return RadicalInverseSpecialized<4523>(a);
    case 615:
        return RadicalInverseSpecialized<4547>(a);
    case 616:
        return RadicalInverseSpecialized<4549>(a);
    case 617:
        return RadicalInverseSpecialized<4561>(a);
    case 618:
        return RadicalInverseSpecialized<4567>(a);
    case 619:
        return RadicalInverseSpecialized<4583>(a);
    case 620:
        return RadicalInverseSpecialized<4591>(a);
    case 621:
        return RadicalInverseSpecialized<4597>(a);
    case 622:
        return RadicalInverseSpecialized<4603>(a);
    case 623:
        return RadicalInverseSpecialized<4621>(a);
    case 624:
        return RadicalInverseSpecialized<4637>(a);
    case 625:
        return RadicalInverseSpecialized<4639>(a);
    case 626:
        return RadicalInverseSpecialized<4643>(a);
    case 627:
        return RadicalInverseSpecialized<4649>(a);
    case 628:
        return RadicalInverseSpecialized<4651>(a);
    case 629:
        return RadicalInverseSpecialized<4657>(a);
    case 630:
        return RadicalInverseSpecialized<4663>(a);
    case 631:
        return RadicalInverseSpecialized<4673>(a);
    case 632:
        return RadicalInverseSpecialized<4679>(a);
    case 633:
        return RadicalInverseSpecialized<4691>(a);
    case 634:
        return RadicalInverseSpecialized<4703>(a);
    case 635:
        return RadicalInverseSpecialized<4721>(a);
    case 636:
        return RadicalInverseSpecialized<4723>(a);
    case 637:
        return RadicalInverseSpecialized<4729>(a);
    case 638:
        return RadicalInverseSpecialized<4733>(a);
    case 639:
        return RadicalInverseSpecialized<4751>(a);
    case 640:
        return RadicalInverseSpecialized<4759>(a);
    case 641:
        return RadicalInverseSpecialized<4783>(a);
    case 642:
        return RadicalInverseSpecialized<4787>(a);
    case 643:
        return RadicalInverseSpecialized<4789>(a);
    case 644:
        return RadicalInverseSpecialized<4793>(a);
    case 645:
        return RadicalInverseSpecialized<4799>(a);
    case 646:
        return RadicalInverseSpecialized<4801>(a);
    case 647:
        return RadicalInverseSpecialized<4813>(a);
    case 648:
        return RadicalInverseSpecialized<4817>(a);
    case 649:
        return RadicalInverseSpecialized<4831>(a);
    case 650:
        return RadicalInverseSpecialized<4861>(a);
    case 651:
        return RadicalInverseSpecialized<4871>(a);
    case 652:
        return RadicalInverseSpecialized<4877>(a);
    case 653:
        return RadicalInverseSpecialized<4889>(a);
    case 654:
        return RadicalInverseSpecialized<4903>(a);
    case 655:
        return RadicalInverseSpecialized<4909>(a);
    case 656:
        return RadicalInverseSpecialized<4919>(a);
    case 657:
        return RadicalInverseSpecialized<4931>(a);
    case 658:
        return RadicalInverseSpecialized<4933>(a);
    case 659:
        return RadicalInverseSpecialized<4937>(a);
    case 660:
        return RadicalInverseSpecialized<4943>(a);
    case 661:
        return RadicalInverseSpecialized<4951>(a);
    case 662:
        return RadicalInverseSpecialized<4957>(a);
    case 663:
        return RadicalInverseSpecialized<4967>(a);
    case 664:
        return RadicalInverseSpecialized<4969>(a);
    case 665:
        return RadicalInverseSpecialized<4973>(a);
    case 666:
        return RadicalInverseSpecialized<4987>(a);
    case 667:
        return RadicalInverseSpecialized<4993>(a);
    case 668:
        return RadicalInverseSpecialized<4999>(a);
    case 669:
        return RadicalInverseSpecialized<5003>(a);
    case 670:
        return RadicalInverseSpecialized<5009>(a);
    case 671:
        return RadicalInverseSpecialized<5011>(a);
    case 672:
        return RadicalInverseSpecialized<5021>(a);
    case 673:
        return RadicalInverseSpecialized<5023>(a);
    case 674:
        return RadicalInverseSpecialized<5039>(a);
    case 675:
        return RadicalInverseSpecialized<5051>(a);
    case 676:
        return RadicalInverseSpecialized<5059>(a);
    case 677:
        return RadicalInverseSpecialized<5077>(a);
    case 678:
        return RadicalInverseSpecialized<5081>(a);
    case 679:
        return RadicalInverseSpecialized<5087>(a);
    case 680:
        return RadicalInverseSpecialized<5099>(a);
    case 681:
        return RadicalInverseSpecialized<5101>(a);
    case 682:
        return RadicalInverseSpecialized<5107>(a);
    case 683:
        return RadicalInverseSpecialized<5113>(a);
    case 684:
        return RadicalInverseSpecialized<5119>(a);
    case 685:
        return RadicalInverseSpecialized<5147>(a);
    case 686:
        return RadicalInverseSpecialized<5153>(a);
    case 687:
        return RadicalInverseSpecialized<5167>(a);
    case 688:
        return RadicalInverseSpecialized<5171>(a);
    case 689:
        return RadicalInverseSpecialized<5179>(a);
    case 690:
        return RadicalInverseSpecialized<5189>(a);
    case 691:
        return RadicalInverseSpecialized<5197>(a);
    case 692:
        return RadicalInverseSpecialized<5209>(a);
    case 693:
        return RadicalInverseSpecialized<5227>(a);
    case 694:
        return RadicalInverseSpecialized<5231>(a);
    case 695:
        return RadicalInverseSpecialized<5233>(a);
    case 696:
        return RadicalInverseSpecialized<5237>(a);
    case 697:
        return RadicalInverseSpecialized<5261>(a);
    case 698:
        return RadicalInverseSpecialized<5273>(a);
    case 699:
        return RadicalInverseSpecialized<5279>(a);
    case 700:
        return RadicalInverseSpecialized<5281>(a);
    case 701:
        return RadicalInverseSpecialized<5297>(a);
    case 702:
        return RadicalInverseSpecialized<5303>(a);
    case 703:
        return RadicalInverseSpecialized<5309>(a);
    case 704:
        return RadicalInverseSpecialized<5323>(a);
    case 705:
        return RadicalInverseSpecialized<5333>(a);
    case 706:
        return RadicalInverseSpecialized<5347>(a);
    case 707:
        return RadicalInverseSpecialized<5351>(a);
    case 708:
        return RadicalInverseSpecialized<5381>(a);
    case 709:
        return RadicalInverseSpecialized<5387>(a);
    case 710:
        return RadicalInverseSpecialized<5393>(a);
    case 711:
        return RadicalInverseSpecialized<5399>(a);
    case 712:
        return RadicalInverseSpecialized<5407>(a);
    case 713:
        return RadicalInverseSpecialized<5413>(a);
    case 714:
        return RadicalInverseSpecialized<5417>(a);
    case 715:
        return RadicalInverseSpecialized<5419>(a);
    case 716:
        return RadicalInverseSpecialized<5431>(a);
    case 717:
        return RadicalInverseSpecialized<5437>(a);
    case 718:
        return RadicalInverseSpecialized<5441>(a);
    case 719:
        return RadicalInverseSpecialized<5443>(a);
    case 720:
        return RadicalInverseSpecialized<5449>(a);
    case 721:
        return RadicalInverseSpecialized<5471>(a);
    case 722:
        return RadicalInverseSpecialized<5477>(a);
    case 723:
        return RadicalInverseSpecialized<5479>(a);
    case 724:
        return RadicalInverseSpecialized<5483>(a);
    case 725:
        return RadicalInverseSpecialized<5501>(a);
    case 726:
        return RadicalInverseSpecialized<5503>(a);
    case 727:
        return RadicalInverseSpecialized<5507>(a);
    case 728:
        return RadicalInverseSpecialized<5519>(a);
    case 729:
        return RadicalInverseSpecialized<5521>(a);
    case 730:
        return RadicalInverseSpecialized<5527>(a);
    case 731:
        return RadicalInverseSpecialized<5531>(a);
    case 732:
        return RadicalInverseSpecialized<5557>(a);
    case 733:
        return RadicalInverseSpecialized<5563>(a);
    case 734:
        return RadicalInverseSpecialized<5569>(a);
    case 735:
        return RadicalInverseSpecialized<5573>(a);
    case 736:
        return RadicalInverseSpecialized<5581>(a);
    case 737:
        return RadicalInverseSpecialized<5591>(a);
    case 738:
        return RadicalInverseSpecialized<5623>(a);
    case 739:
        return RadicalInverseSpecialized<5639>(a);
    case 740:
        return RadicalInverseSpecialized<5641>(a);
    case 741:
        return RadicalInverseSpecialized<5647>(a);
    case 742:
        return RadicalInverseSpecialized<5651>(a);
    case 743:
        return RadicalInverseSpecialized<5653>(a);
    case 744:
        return RadicalInverseSpecialized<5657>(a);
    case 745:
        return RadicalInverseSpecialized<5659>(a);
    case 746:
        return RadicalInverseSpecialized<5669>(a);
    case 747:
        return RadicalInverseSpecialized<5683>(a);
    case 748:
        return RadicalInverseSpecialized<5689>(a);
    case 749:
        return RadicalInverseSpecialized<5693>(a);
    case 750:
        return RadicalInverseSpecialized<5701>(a);
    case 751:
        return RadicalInverseSpecialized<5711>(a);
    case 752:
        return RadicalInverseSpecialized<5717>(a);
    case 753:
        return RadicalInverseSpecialized<5737>(a);
    case 754:
        return RadicalInverseSpecialized<5741>(a);
    case 755:
        return RadicalInverseSpecialized<5743>(a);
    case 756:
        return RadicalInverseSpecialized<5749>(a);
    case 757:
        return RadicalInverseSpecialized<5779>(a);
    case 758:
        return RadicalInverseSpecialized<5783>(a);
    case 759:
        return RadicalInverseSpecialized<5791>(a);
    case 760:
        return RadicalInverseSpecialized<5801>(a);
    case 761:
        return RadicalInverseSpecialized<5807>(a);
    case 762:
        return RadicalInverseSpecialized<5813>(a);
    case 763:
        return RadicalInverseSpecialized<5821>(a);
    case 764:
        return RadicalInverseSpecialized<5827>(a);
    case 765:
        return RadicalInverseSpecialized<5839>(a);
    case 766:
        return RadicalInverseSpecialized<5843>(a);
    case 767:
        return RadicalInverseSpecialized<5849>(a);
    case 768:
        return RadicalInverseSpecialized<5851>(a);
    case 769:
        return RadicalInverseSpecialized<5857>(a);
    case 770:
        return RadicalInverseSpecialized<5861>(a);
    case 771:
        return RadicalInverseSpecialized<5867>(a);
    case 772:
        return RadicalInverseSpecialized<5869>(a);
    case 773:
        return RadicalInverseSpecialized<5879>(a);
    case 774:
        return RadicalInverseSpecialized<5881>(a);
    case 775:
        return RadicalInverseSpecialized<5897>(a);
    case 776:
        return RadicalInverseSpecialized<5903>(a);
    case 777:
        return RadicalInverseSpecialized<5923>(a);
    case 778:
        return RadicalInverseSpecialized<5927>(a);
    case 779:
        return RadicalInverseSpecialized<5939>(a);
    case 780:
        return RadicalInverseSpecialized<5953>(a);
    case 781:
        return RadicalInverseSpecialized<5981>(a);
    case 782:
        return RadicalInverseSpecialized<5987>(a);
    case 783:
        return RadicalInverseSpecialized<6007>(a);
    case 784:
        return RadicalInverseSpecialized<6011>(a);
    case 785:
        return RadicalInverseSpecialized<6029>(a);
    case 786:
        return RadicalInverseSpecialized<6037>(a);
    case 787:
        return RadicalInverseSpecialized<6043>(a);
    case 788:
        return RadicalInverseSpecialized<6047>(a);
    case 789:
        return RadicalInverseSpecialized<6053>(a);
    case 790:
        return RadicalInverseSpecialized<6067>(a);
    case 791:
        return RadicalInverseSpecialized<6073>(a);
    case 792:
        return RadicalInverseSpecialized<6079>(a);
    case 793:
        return RadicalInverseSpecialized<6089>(a);
    case 794:
        return RadicalInverseSpecialized<6091>(a);
    case 795:
        return RadicalInverseSpecialized<6101>(a);
    case 796:
        return RadicalInverseSpecialized<6113>(a);
    case 797:
        return RadicalInverseSpecialized<6121>(a);
    case 798:
        return RadicalInverseSpecialized<6131>(a);
    case 799:
        return RadicalInverseSpecialized<6133>(a);
    case 800:
        return RadicalInverseSpecialized<6143>(a);
    case 801:
        return RadicalInverseSpecialized<6151>(a);
    case 802:
        return RadicalInverseSpecialized<6163>(a);
    case 803:
        return RadicalInverseSpecialized<6173>(a);
    case 804:
        return RadicalInverseSpecialized<6197>(a);
    case 805:
        return RadicalInverseSpecialized<6199>(a);
    case 806:
        return RadicalInverseSpecialized<6203>(a);
    case 807:
        return RadicalInverseSpecialized<6211>(a);
    case 808:
        return RadicalInverseSpecialized<6217>(a);
    case 809:
        return RadicalInverseSpecialized<6221>(a);
    case 810:
        return RadicalInverseSpecialized<6229>(a);
    case 811:
        return RadicalInverseSpecialized<6247>(a);
    case 812:
        return RadicalInverseSpecialized<6257>(a);
    case 813:
        return RadicalInverseSpecialized<6263>(a);
    case 814:
        return RadicalInverseSpecialized<6269>(a);
    case 815:
        return RadicalInverseSpecialized<6271>(a);
    case 816:
        return RadicalInverseSpecialized<6277>(a);
    case 817:
        return RadicalInverseSpecialized<6287>(a);
    case 818:
        return RadicalInverseSpecialized<6299>(a);
    case 819:
        return RadicalInverseSpecialized<6301>(a);
    case 820:
        return RadicalInverseSpecialized<6311>(a);
    case 821:
        return RadicalInverseSpecialized<6317>(a);
    case 822:
        return RadicalInverseSpecialized<6323>(a);
    case 823:
        return RadicalInverseSpecialized<6329>(a);
    case 824:
        return RadicalInverseSpecialized<6337>(a);
    case 825:
        return RadicalInverseSpecialized<6343>(a);
    case 826:
        return RadicalInverseSpecialized<6353>(a);
    case 827:
        return RadicalInverseSpecialized<6359>(a);
    case 828:
        return RadicalInverseSpecialized<6361>(a);
    case 829:
        return RadicalInverseSpecialized<6367>(a);
    case 830:
        return RadicalInverseSpecialized<6373>(a);
    case 831:
        return RadicalInverseSpecialized<6379>(a);
    case 832:
        return RadicalInverseSpecialized<6389>(a);
    case 833:
        return RadicalInverseSpecialized<6397>(a);
    case 834:
        return RadicalInverseSpecialized<6421>(a);
    case 835:
        return RadicalInverseSpecialized<6427>(a);
    case 836:
        return RadicalInverseSpecialized<6449>(a);
    case 837:
        return RadicalInverseSpecialized<6451>(a);
    case 838:
        return RadicalInverseSpecialized<6469>(a);
    case 839:
        return RadicalInverseSpecialized<6473>(a);
    case 840:
        return RadicalInverseSpecialized<6481>(a);
    case 841:
        return RadicalInverseSpecialized<6491>(a);
    case 842:
        return RadicalInverseSpecialized<6521>(a);
    case 843:
        return RadicalInverseSpecialized<6529>(a);
    case 844:
        return RadicalInverseSpecialized<6547>(a);
    case 845:
        return RadicalInverseSpecialized<6551>(a);
    case 846:
        return RadicalInverseSpecialized<6553>(a);
    case 847:
        return RadicalInverseSpecialized<6563>(a);
    case 848:
        return RadicalInverseSpecialized<6569>(a);
    case 849:
        return RadicalInverseSpecialized<6571>(a);
    case 850:
        return RadicalInverseSpecialized<6577>(a);
    case 851:
        return RadicalInverseSpecialized<6581>(a);
    case 852:
        return RadicalInverseSpecialized<6599>(a);
    case 853:
        return RadicalInverseSpecialized<6607>(a);
    case 854:
        return RadicalInverseSpecialized<6619>(a);
    case 855:
        return RadicalInverseSpecialized<6637>(a);
    case 856:
        return RadicalInverseSpecialized<6653>(a);
    case 857:
        return RadicalInverseSpecialized<6659>(a);
    case 858:
        return RadicalInverseSpecialized<6661>(a);
    case 859:
        return RadicalInverseSpecialized<6673>(a);
    case 860:
        return RadicalInverseSpecialized<6679>(a);
    case 861:
        return RadicalInverseSpecialized<6689>(a);
    case 862:
        return RadicalInverseSpecialized<6691>(a);
    case 863:
        return RadicalInverseSpecialized<6701>(a);
    case 864:
        return RadicalInverseSpecialized<6703>(a);
    case 865:
        return RadicalInverseSpecialized<6709>(a);
    case 866:
        return RadicalInverseSpecialized<6719>(a);
    case 867:
        return RadicalInverseSpecialized<6733>(a);
    case 868:
        return RadicalInverseSpecialized<6737>(a);
    case 869:
        return RadicalInverseSpecialized<6761>(a);
    case 870:
        return RadicalInverseSpecialized<6763>(a);
    case 871:
        return RadicalInverseSpecialized<6779>(a);
    case 872:
        return RadicalInverseSpecialized<6781>(a);
    case 873:
        return RadicalInverseSpecialized<6791>(a);
    case 874:
        return RadicalInverseSpecialized<6793>(a);
    case 875:
        return RadicalInverseSpecialized<6803>(a);
    case 876:
        return RadicalInverseSpecialized<6823>(a);
    case 877:
        return RadicalInverseSpecialized<6827>(a);
    case 878:
        return RadicalInverseSpecialized<6829>(a);
    case 879:
        return RadicalInverseSpecialized<6833>(a);
    case 880:
        return RadicalInverseSpecialized<6841>(a);
    case 881:
        return RadicalInverseSpecialized<6857>(a);
    case 882:
        return RadicalInverseSpecialized<6863>(a);
    case 883:
        return RadicalInverseSpecialized<6869>(a);
    case 884:
        return RadicalInverseSpecialized<6871>(a);
    case 885:
        return RadicalInverseSpecialized<6883>(a);
    case 886:
        return RadicalInverseSpecialized<6899>(a);
    case 887:
        return RadicalInverseSpecialized<6907>(a);
    case 888:
        return RadicalInverseSpecialized<6911>(a);
    case 889:
        return RadicalInverseSpecialized<6917>(a);
    case 890:
        return RadicalInverseSpecialized<6947>(a);
    case 891:
        return RadicalInverseSpecialized<6949>(a);
    case 892:
        return RadicalInverseSpecialized<6959>(a);
    case 893:
        return RadicalInverseSpecialized<6961>(a);
    case 894:
        return RadicalInverseSpecialized<6967>(a);
    case 895:
        return RadicalInverseSpecialized<6971>(a);
    case 896:
        return RadicalInverseSpecialized<6977>(a);
    case 897:
        return RadicalInverseSpecialized<6983>(a);
    case 898:
        return RadicalInverseSpecialized<6991>(a);
    case 899:
        return RadicalInverseSpecialized<6997>(a);
    case 900:
        return RadicalInverseSpecialized<7001>(a);
    case 901:
        return RadicalInverseSpecialized<7013>(a);
    case 902:
        return RadicalInverseSpecialized<7019>(a);
    case 903:
        return RadicalInverseSpecialized<7027>(a);
    case 904:
        return RadicalInverseSpecialized<7039>(a);
    case 905:
        return RadicalInverseSpecialized<7043>(a);
    case 906:
        return RadicalInverseSpecialized<7057>(a);
    case 907:
        return RadicalInverseSpecialized<7069>(a);
    case 908:
        return RadicalInverseSpecialized<7079>(a);
    case 909:
        return RadicalInverseSpecialized<7103>(a);
    case 910:
        return RadicalInverseSpecialized<7109>(a);
    case 911:
        return RadicalInverseSpecialized<7121>(a);
    case 912:
        return RadicalInverseSpecialized<7127>(a);
    case 913:
        return RadicalInverseSpecialized<7129>(a);
    case 914:
        return RadicalInverseSpecialized<7151>(a);
    case 915:
        return RadicalInverseSpecialized<7159>(a);
    case 916:
        return RadicalInverseSpecialized<7177>(a);
    case 917:
        return RadicalInverseSpecialized<7187>(a);
    case 918:
        return RadicalInverseSpecialized<7193>(a);
    case 919:
        return RadicalInverseSpecialized<7207>(a);
    case 920:
        return RadicalInverseSpecialized<7211>(a);
    case 921:
        return RadicalInverseSpecialized<7213>(a);
    case 922:
        return RadicalInverseSpecialized<7219>(a);
    case 923:
        return RadicalInverseSpecialized<7229>(a);
    case 924:
        return RadicalInverseSpecialized<7237>(a);
    case 925:
        return RadicalInverseSpecialized<7243>(a);
    case 926:
        return RadicalInverseSpecialized<7247>(a);
    case 927:
        return RadicalInverseSpecialized<7253>(a);
    case 928:
        return RadicalInverseSpecialized<7283>(a);
    case 929:
        return RadicalInverseSpecialized<7297>(a);
    case 930:
        return RadicalInverseSpecialized<7307>(a);
    case 931:
        return RadicalInverseSpecialized<7309>(a);
    case 932:
        return RadicalInverseSpecialized<7321>(a);
    case 933:
        return RadicalInverseSpecialized<7331>(a);
    case 934:
        return RadicalInverseSpecialized<7333>(a);
    case 935:
        return RadicalInverseSpecialized<7349>(a);
    case 936:
        return RadicalInverseSpecialized<7351>(a);
    case 937:
        return RadicalInverseSpecialized<7369>(a);
    case 938:
        return RadicalInverseSpecialized<7393>(a);
    case 939:
        return RadicalInverseSpecialized<7411>(a);
    case 940:
        return RadicalInverseSpecialized<7417>(a);
    case 941:
        return RadicalInverseSpecialized<7433>(a);
    case 942:
        return RadicalInverseSpecialized<7451>(a);
    case 943:
        return RadicalInverseSpecialized<7457>(a);
    case 944:
        return RadicalInverseSpecialized<7459>(a);
    case 945:
        return RadicalInverseSpecialized<7477>(a);
    case 946:
        return RadicalInverseSpecialized<7481>(a);
    case 947:
        return RadicalInverseSpecialized<7487>(a);
    case 948:
        return RadicalInverseSpecialized<7489>(a);
    case 949:
        return RadicalInverseSpecialized<7499>(a);
    case 950:
        return RadicalInverseSpecialized<7507>(a);
    case 951:
        return RadicalInverseSpecialized<7517>(a);
    case 952:
        return RadicalInverseSpecialized<7523>(a);
    case 953:
        return RadicalInverseSpecialized<7529>(a);
    case 954:
        return RadicalInverseSpecialized<7537>(a);
    case 955:
        return RadicalInverseSpecialized<7541>(a);
    case 956:
        return RadicalInverseSpecialized<7547>(a);
    case 957:
        return RadicalInverseSpecialized<7549>(a);
    case 958:
        return RadicalInverseSpecialized<7559>(a);
    case 959:
        return RadicalInverseSpecialized<7561>(a);
    case 960:
        return RadicalInverseSpecialized<7573>(a);
    case 961:
        return RadicalInverseSpecialized<7577>(a);
    case 962:
        return RadicalInverseSpecialized<7583>(a);
    case 963:
        return RadicalInverseSpecialized<7589>(a);
    case 964:
        return RadicalInverseSpecialized<7591>(a);
    case 965:
        return RadicalInverseSpecialized<7603>(a);
    case 966:
        return RadicalInverseSpecialized<7607>(a);
    case 967:
        return RadicalInverseSpecialized<7621>(a);
    case 968:
        return RadicalInverseSpecialized<7639>(a);
    case 969:
        return RadicalInverseSpecialized<7643>(a);
    case 970:
        return RadicalInverseSpecialized<7649>(a);
    case 971:
        return RadicalInverseSpecialized<7669>(a);
    case 972:
        return RadicalInverseSpecialized<7673>(a);
    case 973:
        return RadicalInverseSpecialized<7681>(a);
    case 974:
        return RadicalInverseSpecialized<7687>(a);
    case 975:
        return RadicalInverseSpecialized<7691>(a);
    case 976:
        return RadicalInverseSpecialized<7699>(a);
    case 977:
        return RadicalInverseSpecialized<7703>(a);
    case 978:
        return RadicalInverseSpecialized<7717>(a);
    case 979:
        return RadicalInverseSpecialized<7723>(a);
    case 980:
        return RadicalInverseSpecialized<7727>(a);
    case 981:
        return RadicalInverseSpecialized<7741>(a);
    case 982:
        return RadicalInverseSpecialized<7753>(a);
    case 983:
        return RadicalInverseSpecialized<7757>(a);
    case 984:
        return RadicalInverseSpecialized<7759>(a);
    case 985:
        return RadicalInverseSpecialized<7789>(a);
    case 986:
        return RadicalInverseSpecialized<7793>(a);
    case 987:
        return RadicalInverseSpecialized<7817>(a);
    case 988:
        return RadicalInverseSpecialized<7823>(a);
    case 989:
        return RadicalInverseSpecialized<7829>(a);
    case 990:
        return RadicalInverseSpecialized<7841>(a);
    case 991:
        return RadicalInverseSpecialized<7853>(a);
    case 992:
        return RadicalInverseSpecialized<7867>(a);
    case 993:
        return RadicalInverseSpecialized<7873>(a);
    case 994:
        return RadicalInverseSpecialized<7877>(a);
    case 995:
        return RadicalInverseSpecialized<7879>(a);
    case 996:
        return RadicalInverseSpecialized<7883>(a);
    case 997:
        return RadicalInverseSpecialized<7901>(a);
    case 998:
        return RadicalInverseSpecialized<7907>(a);
    case 999:
        return RadicalInverseSpecialized<7919>(a);
    case 1000:
        return RadicalInverseSpecialized<7927>(a);
    case 1001:
        return RadicalInverseSpecialized<7933>(a);
    case 1002:
        return RadicalInverseSpecialized<7937>(a);
    case 1003:
        return RadicalInverseSpecialized<7949>(a);
    case 1004:
        return RadicalInverseSpecialized<7951>(a);
    case 1005:
        return RadicalInverseSpecialized<7963>(a);
    case 1006:
        return RadicalInverseSpecialized<7993>(a);
    case 1007:
        return RadicalInverseSpecialized<8009>(a);
    case 1008:
        return RadicalInverseSpecialized<8011>(a);
    case 1009:
        return RadicalInverseSpecialized<8017>(a);
    case 1010:
        return RadicalInverseSpecialized<8039>(a);
    case 1011:
        return RadicalInverseSpecialized<8053>(a);
    case 1012:
        return RadicalInverseSpecialized<8059>(a);
    case 1013:
        return RadicalInverseSpecialized<8069>(a);
    case 1014:
        return RadicalInverseSpecialized<8081>(a);
    case 1015:
        return RadicalInverseSpecialized<8087>(a);
    case 1016:
        return RadicalInverseSpecialized<8089>(a);
    case 1017:
        return RadicalInverseSpecialized<8093>(a);
    case 1018:
        return RadicalInverseSpecialized<8101>(a);
    case 1019:
        return RadicalInverseSpecialized<8111>(a);
    case 1020:
        return RadicalInverseSpecialized<8117>(a);
    case 1021:
        return RadicalInverseSpecialized<8123>(a);
    case 1022:
        return RadicalInverseSpecialized<8147>(a);
    case 1023:
        return RadicalInverseSpecialized<8161>(a);
    default:
        return 0;
    }
}