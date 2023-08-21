#ifndef CPBRT_CORE_MICROFACET_H
#define CPBRT_CORE_MICROFACET_H

#include "geometry.h"

// A microfacet distribution describes the microgeometry of the surface in terms of the slope and
// azimuthal orientation of V-shaped microfacets. 
class MicrofacetDistribution {
protected:
    // TODO.
    const bool sampleVisibleArea;

    MicrofacetDistribution(bool sampleVisibleArea) : sampleVisibleArea(sampleVisibleArea) {}

public:

    // Microfacet distribution function that gives the total differential area
    // of microfacets that have the given normal (in the reflection coordinate system).
    virtual Float D(const Vector3f &wh) const = 0;

    // Smith's geometric attenuation factor G(wo, wi) gives the fraction of microfacets in a
    // differential area dA that are visible from both directions wo and wi. We assume that
    // visibility is more likely the higher up a given point on a microfacet is.
    Float G(const Vector3f &wo, const Vector3f &wi) const {
        return 1 / (1 + Lambda(wo) + Lambda(wi));
    }

    // 0 <= G1(w) <= 1 is Smith's masking-shadowing function that gives the fraction of
    // normalized and projected microfacet area that is visible from the direction w. If
    // we let A+(w) denote the projected area of forward facing microfacets in the direction w
    // A-(w) the projected area of backfacing microfacets, then G1(w) = [A+(w) - A-(w)]/A+(w)
    // is the ratio of visible forward-facing microfacet area to total forward-facing microfacet
    // area. 
    Float G1(const Vector3f &w) const {
        return 1 / (1 + Lambda(w));
    }

    // Lambda(w) is the ratio of masked forward-facing microfacet area per visible forward-facing
    // microfacet area. Depends on the specific distribution, so no base-class implementation exists.
    virtual Float Lambda(const Vector3f &w) const = 0;

    // TODO.
    virtual Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const = 0;

    // TODO.
    Float Pdf(const Vector3f &wo, const Vector3f &wh) const;
};

// Beckmann-Spizzichino microfacet distribution.
class BeckmannDistribution : public MicrofacetDistribution {
private:
    // Parameters of the distribution. Alpha denotes the root mean square (RMS) slope (polar angle)
    // of microfacets with a given normal. alphaX is the RMS slope corresponding to microfacets
    // whose azimuthal orientation is perpendicular to the x-axis of the reflection coordinate
    // system; alphaY is the one for microfacets whose azimuthal orientation is perpendicular to the
    // y-axis. 
    const Float alphaX, alphaY;
public:
    BeckmannDistribution(Float alphaX, Float alphaY, bool sampleVisibleArea = true)
        : MicrofacetDistribution(sampleVisibleArea), 
          alphaX(std::max(Float(0.001), alphaX)),
          alphaY(std::max(Float(0.001), alphaY))
    {}

    // Beckmann-Spizzichino microfacet distribution function. Anisotropic in general (dependent on
    // azimuthal angle phi), isotropic in the case where alphaX = alphaY.
    Float D(const Vector3f &wh) const;

    Float Lambda(const Vector3f &w) const;

    // TODO.
    Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;

    // TODO.
    Float Pdf(const Vector3f &wo, const Vector3f &wh) const;

    // Converts a [0,1] roughness value, where 0 is perfect specular and 1 is very rough,
    // to root mean square slope, alpha.
    static Float RoughnessToAlpha(Float roughness) {
        roughness = std::max(roughness, (Float) 1e-3);
        Float x = std::log(roughness);
        return 1.62142f + 0.819955f*x + 0.1734f*x*x + 0.0171201f*x*x*x + 0.000640711f*x*x*x*x;
    }
};

// Trowbridge-Reitz microfacet distribution.
class TrowbridgeReitzDistribution : public MicrofacetDistribution {
private:
    // Parameters of the distribution. Alpha X and Y control the roughness of the
    // surface. RoughnessToAlpha() maps a roughness value in the typical [0,1] range
    // to alpha X and Y values. 
    const Float alphaX, alphaY;

public:
    TrowbridgeReitzDistribution(Float alphaX, Float alphaY, bool sampleVisibleArea = true)
     : MicrofacetDistribution(sampleVisibleArea),
       alphaX(std::max(Float(0.001), alphaX)),
       alphaY(std::max(Float(0.001), alphaY))
    {}

    // Trowbridge-Reitz microfacet distribution function. Anisotropic in general (dependent on
    // azimuthal angle phi), isotropic in the case where alphaX = alphaY.
    Float D(const Vector3f &wh) const;

    Vector3f Sample_wh(const Vector3f &wo, const Point2f &u) const;

    static inline Float RoughnessToAlpha(Float roughness);

private:
    Float Lambda(const Vector3f &w) const;
};

// Maps a roughness value in the range [0,1] to a value for one of the Trowbridge-Reitz
// alpha parameters.
inline Float TrowbridgeReitzDistribution::RoughnessToAlpha(Float roughness) {
    // Not explained in PBRT.
    roughness = std::max(roughness, (Float) 1e-3);
    Float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}

// Disney's microfacet distribution.
class DisneyMicrofacetDistribution : public TrowbridgeReitzDistribution {
public:
    DisneyMicrofacetDistribution(Float alphaX, Float alphaY)
     : TrowbridgeReitzDistribution(alphaX, alphaY)
    {}

    // Smith's geometric attenuation factor G(wo, wi) gives the fraction of microfacets in a
    // differential area dA that are visible from both directions wo and wi. We assume that
    // visibility is more likely the higher up a given point on a microfacet is.
    Float G(const Vector3f &wo, const Vector3f &wi) const {
        // 0 <= G1(w) <= 1 is Smith's masking-shadowing function that gives the fraction of
        // normalized and projected microfacet area that is visible from the direction w.
        //
        // Separable: G(wo, wi) = G1(wo)G1(wi)
        //            G(wo, wi)/G1(wo) = G1(wi)
        //            G(wo, wi)/G1(wi) = G1(wo)
        return G1(wo) * G1(wi);
    }
};

#endif // CPBRT_CORE_MICROFACET_H