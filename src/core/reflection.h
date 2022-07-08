#ifndef CPBRT_CORE_REFLECTION_H
#define CPBRT_CORE_REFLECTION_H

#include "cpbrt.h"
#include "geometry.h"
#include "shape.h"
#include "spectrum.h"
#include "microfacet.h"

// The shading coordinate system is defined by the orthonormal basis {s, t, n} = {x, y, z},
// where s and t are 2 orthogonal vectors tangent to the shaded point and n is the normal
// of the surface at this point.
//
// This coordinate system defines the spherical coordinates of a vector w as (theta, phi),
// where theta is measured from the z axis and phi is measured about the z axis and from the
// x axis to the projection of w onto the xy plane (the length of the projection is sin(theta)).
//
// Vector w is expected to be expressed as a linear combination of {s, t, n}, normalized, and
// outward facing (even if it is an incident direction). When it is an incident direction, it
// will always be in the same hemisphere as the normal n. 

inline Float CosTheta(const Vector3f &w) {
    return w.z;
}

inline Float Cos2Theta(const Vector3f &w) {
    return w.z * w.z;
}

inline Float AbsCosTheta(const Vector3f &w) {
    return std::abs(w.z);
}

inline Float Sin2Theta(const Vector3f &w) {
    // Pythagorean identity.
    return std::max((Float) 0, (Float) 1 - Cos2Theta(w));
}

inline Float SinTheta(const Vector3f &w) {
    return std::sqrt(Sin2Theta(w));
}

inline Float TanTheta(const Vector3f &w) {
    // Trigonometric identity.
    return SinTheta(w) / CosTheta(w);
}

inline Float Tan2Theta(const Vector3f &w) {
    // Trigonometric identity.
    return Sin2Theta(w) / Cos2Theta(w);
}

inline Float CosPhi(const Vector3f &w) {
    // The length of the projection of w onto the xy plane where phi is measured is given
    // by sin(theta).
    Float sinTheta = SinTheta(w);

    // Trigonometric identity (the projection of w is the hypotenuse of the triangle).
    return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}

inline Float Cos2Phi(const Vector3f &w) {
    return CosPhi(w) * CosPhi(w);
}

inline Float SinPhi(const Vector3f &w) {
    // The length of the projection of w onto the xy plane where phi is measured is given
    // by sin(theta).
    Float sinTheta = SinTheta(w);

    // Trigonometric identity (the projection of w is the hypotenuse of the triangle).
    return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}

inline Float Sin2Phi(const Vector3f &w) {
    return SinPhi(w) * SinPhi(w);
}

// Computes the cosine of the angle between vectors wa and wb.
inline Float CosDPhi(const Vector3f &wa, const Vector3f &wb) {
    // When wa and wb aren't normalized, their dot product gives ||wa||*||wb||*cos(alpha),
    // where alpha is the angle between them. So it must be divided by ||wa||*||wb||.
    //
    // (The dot product of wa with itself gives the square of its norm ||wa||.) 
    return Clamp(
        (wa.x*wb.x + wa.y*wb.y) / std::sqrt((wa.x*wa.x + wa.y*wa.y) * (wb.x*wb.x + wb.y*wb.y)),
        -1, 1
    );
}

inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) {
    return w.z * wp.z > 0;
}

inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) {
    return -wo + 2*Dot(wo, n)* n;
}

inline bool Refract(const Vector3f &wi, const Normal3f &n, Float eta, Vector3f *wt) {
    // ThetaI is the angle of incidence, measured with respect to the normal in reflection
    // space.
    Float cosThetaI = Dot(n, wi);
    Float sin2ThetaI = std::max((Float) 0.f, 1.f - cosThetaI * cosThetaI);
    
    // ThetaT is the refraction/transmission angle determined by Snell's law.
    Float sin2ThetaT = eta * eta * sin2ThetaI;

    if (sin2ThetaT >= 1) {
        // Total internal reflection. When light travels from a medium to another with
        // lower refraction index and it does so past a critical angle, tending to
        // graze the boundary, Snell's law doesn't have a solution and refraction can't
        // occur, so light gets reflected instead back into the incident medium.
        return false;
    }
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);

    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
    return true;
}

// Evaluates the Fresnel reflectance equation between 2 dielectric media, assuming that light is
// unpolarized. cosThetaI is the angle of incidence measured from the normal; etaI is the refraction
// index of the medium that light is traveling through before reaching the interface with the
// other medium, whose refraction index etaT is. (A refraction index is a property of the medium:
// the ratio of the speed of light in a vacuum to the speed of light through the medium. Refraction
// indices for dielectrics are assumed to be real numbers.)
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT);

// Evaluates the general Fresnel equation, but assuming that the transmission medium is a dielectric.
// Refraction indices are spectra because they are wavelength-dependent. cosThetaI is the angle of
// incidence measured from the normal; etaI and etaT are the indices of refraction of the incident
// and transmission media, respectively; and k is the imaginary part of the index of refraction of
// the incident medium (assumed to be a conductor) also known as absorption coefficient.
Spectrum FrConductor(Float cosThetaI, const Spectrum &etaI, const Spectrum &etaT, const Spectrum &k);

class Fresnel {
public:
    // Evaluates the Fresnel reflectance equation for the input incidence angle (measured from the
    // normal).
    virtual Spectrum Evaluate(Float cosI) const = 0;
};

class FresnelDielectric : public Fresnel {
private:
    // Index of refraction of the incident dielectric medium.
    Float etaI;

    // Index of refraction of the transmission dielectric medium.
    Float etaT;

public:
    FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) {}

    Spectrum Evaluate(Float cosThetaI) const;
};

class FresnelConductor : public Fresnel {
private:
    // Index of refraction of the incident medium.
    Spectrum etaI;

    // Index of refraction of the transmission medium (assumed to be a dielectric).
    Spectrum etaT;

    // Imaginary part of the index of refraction of the incident medium, also known as absorption
    // coefficient.
    Spectrum k;

public:
    FresnelConductor(const Spectrum &etaI, const Spectrum &etaT, const Spectrum &k)
        : etaI(etaI), etaT(etaT), k(k)
    {}

    Spectrum Evaluate(Float cosThetaI) const;
};

// A Fresnel interface that reflects all the incident light in its entirety: no absorption, no
// transmission.
class FresnelNoOp : public Fresnel {
public:
    Spectrum Evaluate(Float) const {
        return Spectrum(1.);
    }
};

enum BxDFType {
    BSDF_REFLECTION   = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE      = 1 << 2,
    BSDF_GLOSSY       = 1 << 3,
    BSDF_SPECULAR     = 1 << 4,
    BSDF_ALL          = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR
};

// Interface for BRDF and BTDF.
class BxDF {
public:
    const BxDFType type;

    BxDF(BxDFType type) : type(type) {}

    bool MatchesFlags(BxDFType t) const {
        return (type & t) == type;
    }

    // Computes the spectral distribution of a radiometric quantity over wavelength in the
    // given outgoing direction that corresponds to the given incident direction.
    virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;

    // Computes the spectral distribution of a radiometric quantity over wavelength in the
    // given outgoing direction, computing also the (possibly unique) corresponding incident
    // direction. Vectors wo and wi are expressed with respect to the local coordinate system.
    virtual Spectrum Sample_f(
        const Vector3f &wo,
        Vector3f *wi,
        const Point2f &sample,
        Float *pdf,
        BxDFType *sampledType = nullptr
    ) const;

    // Computes the hemispherical-directional reflectance of the surface in the given outgoing
    // direction considering the entire hemisphere of incident directions. 
    //
    // Reflectance is the spectral distribution that describes the reflective (scattering)
    // behavior of the surface. Reflectance alone does not compute radiance, but it dictates
    // the spectral distribution of it.
    //
    // nSamples and samples are only used when the implementation can't compute it exactly using
    // a closed form and must approximate it instead.
    virtual Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples) const;

    // Computes the hemispherical-hemispherical reflectance of the surface.
    // TODO: explain.
    //
    // nSamples and samples are only used when the implementation can't compute it exactly using
    // a closed form and must approximate it instead.
    virtual Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const;

    virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
};

class ScaledBxDF : public BxDF {
public:
    ScaledBxDF(BxDF *bxdf, const Spectrum &scale) 
        : BxDF(BxDFType(bxdf->type)), bxdf(bxdf), scale(scale)
    {}

    // Scales the spectral distribution of the wrapped BxDF by scale.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    // Scales the spectral distribution of the wrapped BxDF by scale, computing also the (possibly
    // unique) corresponding incident direction.
    Spectrum Sample_f(
        const Vector3f &wo,
        Vector3f *wi,
        const Point2f &sample,
        Float *pdf,
        BxDFType *sampledType = nullptr
    ) const;

    // Scales the hemispherical-directional reflectance of the wrapped BxDF by scale.
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples) const;

    // Scales the hemispherical-hemispherical reflectance of the wrapped BxDF by scale.
    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const;
  
private:
    BxDF *bxdf;
    Spectrum scale;
};

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDF : public BxDF {
public:
    LayeredBxDF() = default;

    LayeredBxDF(
        TopBxDF top,
        BottomBxDF bottom,
        Float thickness,
        const SampledSpectrum &albedo,
        Float g,
        int maxDepth,
        int nSamples
    ) : BxDF(LayeredBxDF<TopBxDF, BottomBxDF>::DetermineType(top, bottom, albedo)),
        top(top),
        bottom(bottom),
        thickness(std::max(thickness, std::numeric_limits<Float>::min())),
        albedo(albedo),
        g(g),
        maxDepth(maxDepth),
        nSamples(nSamples)
    {}

    // Determines the type (flags) of the layered BxDF from the types of the top and bottom BxDFs.
    static BxDFType DetermineType(TopBxDF top, BottomBxDF bottom, bool hasAlbedo) {
        BxDFType topType = top.type;
        BxDFType bottomType = bottom.type;

        BxDFType type = BxDFType::BSDF_REFLECTION;

        if(topType & BxDFType::BSDF_SPECULAR) {
            type |= BxDFType::BSDF_SPECULAR;
        }

        // Either one can be diffuse (or glossy) for the layered one to be diffuse (glossy).
        if (hasAlbedo || topType & BxDFType::BSDF_DIFFUSE || bottomType & BxDFType::BSDF_DIFFUSE) {
            type |= BxDFType::BSDF_DIFFUSE;
        } else if (topType & BxDFType::BSDF_GLOSSY || bottomType & BxDFType::BSDF_GLOSSY) {
            type |= BxDFType::BSDF_GLOSSY;
        }

        // Both need to be transmissive for the layered one to be transmissive.
        if (topType & BxDFType::BSDF_TRANSMISSION && bottomType & BxDFType::BSDF_TRANSMISSION) {
            type |= BxDFType::BSDF_TRANSMISSION;
        }

        return type;
    }

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

private:
    TopBxDF top;
    BottomBxDF bottom;

    Float thickness;

    Float g;

    // Absorption of the medium. Important for modeling transmission through a possible medium
    // found between the layers.
    SampledSpectrum albedo;

    int maxDepth;
    
    int nSamples;

    // Computes beam transmittance of the medium between layers, if present.
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min()) {
            return 1;
        }
        return FastExp(-std::abs(dz / w.z));
    }
};

template <typename TopBxDF, typename BottomBxDF>
class TopOrBottomBxDF : public BxDF {
public:
    TopOrBottomBxDF() = default;

    TopOrBottomBxDF &operator=(const TopBxDF *t) {
        top = t;
        bottom = nullptr;
        type = t->type;
        return *this;
    }

    TopOrBottomBxDF &operator=(const BottomBxDF *b) {
        bottom = b;
        top = nullptr;
        type = b->type;
        return *this;
    }

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return top ? top->f(wo, wi) : bottom->f(wo, wi);
    }

    Spectrum Sample_f(
        const Vector3f &wo,
        Vector3f *wi,
        const Point2f &sample,
        Float *pdf,
        BxDFType *sampledType = nullptr
    ) const {
        return top ? top->Sample_f(wo, wi, sample, pdf, sampledType) : bottom->Sample_f(wo, wi, sample, pdf, sampledType);
    }

    Float Pdf(const Vector3f &wo, const Vector3f &wi) const {
        return top ? top->Pdf(wo, wi) : bottom->Pdf(wo, wi);
    }

private:
    const TopBxDF *top = nullptr;
    const BottomBxDF *bottom = nullptr;
};

// The Lambertian reflection model models a perfect diffuse surface that reflects incident
// light equally in all directions.
class LambertianReflection : public BxDF {
private:
    // Reflectance. Spectral distribution for all pairs of incident and outgoing directions
    // in the hemisphere. Its value is R=PI*f, where f is a constant BRDF. Since f is a spectral
    // distribution over wavelength, it is referred to as "diffuse color" or "albedo".
    const Spectrum R;

public:
    LambertianReflection(const Spectrum &R)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) 
    {}

    // The reflectance expression includes the BRDF f: R=PI*f.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return R * InvPi;
    }

    // The Lambertian hemispherical-directional reflectance is constant across the entire
    // hemisphere of outgoing directions.
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples) const {
        return R;
    }

    // The Lambertian hemispherical-hemispherical reflectance is constant for all pairs of
    // incident and outgoing directions in the hemisphere.
    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const {
        return R;
    }
};

// TODO: describe.
class LambertianTransmission : public BxDF {
private:
    // Transmittance. Its value is T=PI*f, where f is a constant BTDF.
    const Spectrum T;

public:
    LambertianTransmission(const Spectrum &T)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)), T(T) 
    {}

    // The transmittance expression includes the BTDF f: T=PI*f.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return T * InvPi;
    }

    // TODO: explain.
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples) const {
        return T;
    }

    // TODO: explain.
    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const {
        return T;
    }
};

class SpecularReflection : public BxDF {
private:
    // TODO: ?
    const Spectrum R;

    // Fresnel reflectance equations.
    const Fresnel *fresnel;

public:
    SpecularReflection(const Spectrum &R, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), R(R), fresnel(fresnel)
    {}

    // Computes the spectral distribution of a radiometric quantity over wavelength for an
    // arbitrary pair of outgoing and incident directions. Since there's no chance that an
    // arbitrary pair will satisfy the perfect reflection relation, the returned reflectance
    // is 0.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        // TODO: ?
        return Spectrum(0.f);
    }

    // Samples the BRDF. The only possible incoming direction wi for the input outgoing wo
    // is its reflection about the surface normal. The normal vector doesn't need to be known
    // because it corresponds to the vertical axis in the reflection coordinate system.
    Spectrum Sample_f(
        const Vector3f &wo,
        Vector3f *wi,
        const Point2f &sample,
        Float *pdf,
        BxDFType *sampledType
    ) const;

    // Evaluates the probability that incident direction wi gets sampled for the given outgoing
    // direction wo. Since there's virtually no chance that wi will be the perfect specular
    // reflection direction of wo obtained at random, the probability is 0.
    Float Pdf(const Vector3f &wo, const Vector3f &wi) {
        return 0;
    }
};

class SpecularTransmission : public BxDF {
private:
    // Transmission scale factor.
    // TODO: ?
    const Spectrum T;

    // Incoming medium's refraction index.
    const Float etaA;

    // Transmission medium's refraction index.
    const Float etaB;

    // Fresnel equations.
    const FresnelDielectric fresnel;

    // Did the ray intersecting the point where this BTDF is computed start at the camera
    // or at a light source?
    const TransportMode mode;

public:
    SpecularTransmission(
        const Spectrum &T, Float etaA, Float etaB, TransportMode mode
    ) : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
        T(T),
        etaA(etaA),
        etaB(etaB),
        fresnel(etaA, etaB),
        mode(mode)
    {}

    // BTDF. Computes the spectral distribution of a radiometric quantity over wavelength for an
    // arbitrary pair of transmission and incident directions. Since there's no chance that an
    // arbitrary pair will satisfy the perfect transmission relation, the returned reflectance
    // is 0. A Dirac delta distribution indeed.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f); 
    }

    // Samples the BTDF. The only possible incident direction wi for the input transmission wo
    // is governed by Snell's law of refraction. The normal vector doesn't need to be known
    // because it corresponds to the vertical axis in the reflection coordinate system.
    Spectrum Sample_f(
        const Vector3f &wo, Vector3f *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType
    ) const;
};

// Called just FresnelSpecular in PBRT.
class FresnelSpecularReflectionTransmission : public BxDF {
private:
    // TODO: ?
    const Spectrum R;

    // Transmission scale factor.
    const Spectrum T;

    // Incoming medium's refraction index.
    Float etaA;
    // Transmission medium's refraction index.
    Float etaB;

    // Fresnel equations.
    const FresnelDielectric fresnel;

    // Did the ray intersecting the point where this BTDF is computed start at the camera
    // or at a light source?
    const TransportMode mode;

public:
    FresnelSpecularReflectionTransmission(
        const Spectrum &R, const Spectrum &T, Float etaA, Float etaB, TransportMode mode
    ) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
        R(R), T(T), etaA(etaA), etaB(etaB), fresnel(etaA, etaB), mode(mode)
    { }

    // BSDF. See SpecularReflection::f and SpecularTransmission::f.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const { 
        return Spectrum(0.f);
    }

    // Samples the BSDF, which has a Dirac delta distribution. u is a sample from a uniform distribution
    // that is used as a probability threshold to choose between sampling the BRDF (< u) and the BTDF(> u).
    Spectrum Sample_f(
        const Vector3f &wo, Vector3f *wi, const Point2f &u, Float *pdf, BxDFType *sampledType
    ) const;
};

// Oren-Nayar models diffuse reflection of rough surfaces (which appear brighter as the viewing
// direction approaches the incident direction), which, in general, don't exhibit perfect Lambertian
// reflection (where the brightness of the surface doesn't appear to change as the viewing direction
// changes).
//
// See carlos-lopez-garces.github.io/2021/11/24/oren-nayar-reflectance-model.html for more details.
class OrenNayarReflection : public BxDF {
private:
    // Albedo.
    const Spectrum R;

    // A = 1 - (sigma^2 / 2(sigma^2 + 0.33)), where sigma^2 is the standard deviation of the
    // microfacet angle distribution. With fixed sigma, A is a constant and can be precomputed.
    Float A;

    // B = 0.45sigma^2 / (sigma^2 + 0.09). See A.
    Float B;

public:
    // sigma is the standard deviation of the microfacet angle distribution, in degrees. When sigma=0,
    // all microfacets have the same orientation and the Oren-Nayar reflectance model reduces to the
    // Lambertian model. (Oren-Nayar is, in fact, a generalization of the Lambertian reflectance model.)
    OrenNayarReflection(const Spectrum &R, Float sigma)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {
        
        sigma = Radians(sigma);

        // Variance.
        Float sigma2 = sigma * sigma;

        // Part of the BRDF expression.
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }

    // The Oren-Nayar BRDF is an approximation of the aggregate effect of V-shaped microfacets,
    // where each microfacet exhibits Lambertian reflection.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
};

// Torrance-Sparrow microfacet model. V-shaped cavities of symmetric microfacets with perfectly
// specular reflection. Only the microfacets with a normal equal to the half-angle of a pair of
// directions wo and wi reflect light.
class TorranceSparrowMicrofacetReflection : public BxDF {
private:
    // Single microfacet reflectance.
    const Spectrum R;

    // Distribution of slope and orientation of V-shaped microfacets. The distribution function
    // gives the normalized differential area of microfacets with a given normal wh. Gives also
    // the geometric attenuation factor, GAF, that accounts for masking and shadowing. 
    const MicrofacetDistribution *distribution;

    // Fresnel reflectance. Evaluates the Fresnel reflectance equations that determine the 
    // fraction of light that is reflected (the rest is transmitted or absorbed).
    const Fresnel *fresnel;

public:
    TorranceSparrowMicrofacetReflection(const Spectrum &R, MicrofacetDistribution *distribution, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), R(R), distribution(distribution), fresnel(fresnel)
    {}

    // The Torrance-Sparrow BRDF is given by f(wo, wi) = D(wh)G(wo, wi)Fr(wo)/4cos(thetaO)cos(thetaI),
    // where D(wh) is the microfacet distribution function evaluated at the half-angle that corresponds
    // to wo and wi; G(wo, wi) is the geometric attenuation factor, GAF, that reduces reflected radiance
    // to account for shadowing and masking; and Fr(wo) is the Fresnel reflectance, the fraction of
    // incident radiance that gets reflected.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

    // TODO.
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u, Float *pdf, BxDFType *sampledType) const;
    
    // TODO.
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
};

class BSDF {
private:
    // Geometric normal.
    const Normal3f ng;

    // Orthonormal basis (ss, ts, ns) for shading coordinate system. All are coordinate
    // vectors relative to the world space basis.

    // Shading normal.
    const Normal3f ns;
    // Primary tangent.
    const Vector3f ss;
    // Secondary tangent.
    const Vector3f ts;

    // Component BxDFs at the point.
    int nBxDFs = 0;
    static constexpr int MaxBxDFs = 8;
    BxDF *bxdfs[MaxBxDFs];

    // Private so that it isn't called explicitly using free(). The memory of BSDFs is
    // meant to be managed by a MemoryArena. Allocate BSDFs using ARENA_ALLOC.
    ~BSDF() {}

public:
    // Relative index of refraction at the surface boundary where the point lies. 
    const Float eta;

    BSDF(const SurfaceInteraction &si, Float eta = 1)
    : ns(si.shading.n),
      ng(si.n),
      ss(Normalize(si.shading.dpdu)),
      ts(Cross(ns, ss)),
      eta(eta)
    {}

    // Adds the input BxDF to the collection.
    void Add(BxDF *bxdf) {
        Assert(nBxDFs < MaxBxDFs);
        bxdfs[nBxDFs++] = bxdf;
    }

    int NumComponents(BxDFType flags = BSDF_ALL) const;

    // Change of coordinate from world space to shading space.
    Vector3f WorldToLocal(const Vector3f &v) const {
        // TODO: explain; this is not the change-of-coordinate matrix I'm familiar
        // with or the change of basis.
        return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    }

    // Change of coordinate from shading space to world space.
    Vector3f LocalToWorld(const Vector3f &v) const {
        // TODO: explain; is this the inverse of WorldToLocal computed as the transpose
        // because it's orthogonal?
        return Vector3f(
            ss.x * v.x + ts.x * v.y + ns.x * v.z,
            ss.y * v.x + ts.y * v.y + ns.y * v.z,
            ss.z * v.x + ts.z * v.y + ns.z * v.z
        );
    }

    // Computes the sum of the BRDFs when wo and wi are on the same hemisphere; the
    // sum of the BTDFs otherwise. Only BxDFs of the type specified contribute to the sum.
    // wo and wi are in world space.
    Spectrum f(const Vector3f &woW, const Vector3f &wiW, BxDFType flags = BSDF_ALL) const; 

    // Samples the collection of BxDFs of the input type to compute the spectral distribution
    // of a radiometric quantity over wavelength in the given outgoing direction, computing also
    // the (possibly unique) corresponding incident direction.
    Spectrum Sample_f(
        const Vector3f &woWorld,
        Vector3f *wiWorld,
        const Point2f &u,
        Float *pdf,
        BxDFType type = BSDF_ALL,
        BxDFType *sampledType = nullptr
    ) const;

    // Computes the sum of the hemispherical-directional reflectances of the types specified.
    Spectrum rho(
        const Vector3f &wo,
        int nSamples,
        const Point2f *samples,
        BxDFType flags = BSDF_ALL
    ) const;

    // Computes the sum of the hemispherical-hemispherical reflectances of the types specified.
    Spectrum rho(
        int nSamples,
        const Point2f *samples1,
        const Point2f *samples2,
        BxDFType flags = BSDF_ALL
    ) const;

    Float Pdf(const Vector3f &woWorld, const Vector3f &wiWorld, BxDFType flags = BSDF_ALL) const;
};

inline int BSDF::NumComponents(BxDFType flags) const {
    int num = 0;
    for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++num;
        }
    }
    return num;
}

#endif // CPBRT_CORE_REFLECTION_H