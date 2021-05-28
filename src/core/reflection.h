#ifndef CPBRT_CORE_REFLECTION_H
#define CPBRT_CORE_REFLECTION_H

#include "cpbrt.h"
#include "spectrum.h"

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