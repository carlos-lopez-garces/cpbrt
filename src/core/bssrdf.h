#ifndef CPBRT_CORE_BSSRDF_H
#define CPBRT_CORE_BSSRDF_H

#include "geometry.h"
#include "interaction.h"
#include "reflection.h"

// Computes the 1st moment of the Fresnel reflectance function Fr.
Float FresnelMoment1(Float reciprocalEta);

// Computes the 2nd moment of the Fresnel reflectance function Fr.
Float FresnelMoment2(Float reciprocalEta);

class BSSRDF {
protected:
    const SurfaceInteraction &po;
    Float eta;

public:
    // po is the SurfaceInteraction at the point where exitant radiance is to be computed
    // and eta is the index of refraction of the medium.
    BSSRDF(const SurfaceInteraction &po, Float eta) : po(po), eta(eta) {}

    // Computes the ratio of differential exitant radiance at po in the direction 
    // SurfaceInteraction:wo to the incident differential flux at pi from direction wi.
    virtual Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) = 0;
};

// A separable BSSRDF approximates a BSSRDF as a product of 3 functions:
//
// S(po, wo, pi, wi) ~= (1 - Fr(cos(Theta_o)))Sp(po,pi)Sw(wi)
//
// where Fr is Fresnel transmittance and 1 - Fr(cos(Theta_o)) is the fraction of incident
// radiance that gets transmitted into direction wo when exiting the surface; Sw(wi) is
// the fraction of incident radiance coming from wi and that makes it from pi to po as a 
// result of the influence of the characteristics of the surface boundary at pi; and Sp(po,pi)
// is a spatial distribution that tells how far light travels within the surface.
//
// This approximation effectively decouples the spatial and directional arguments of S.
//
// The separable approximation simplifies the subsurface scattering equation Lo(po,wo) (an
// integral equation) by allowing us to integrate solely with respect to incident direction
// first and then with respect to area only. 
class SeparableBSSRDF : public BSSRDF {
private:
    // Local coordinate frame with origin at po, which corresponds to the shading
    // coordinate system at po.
    const Normal3f ns;
    const Vector3f ss;
    const Vector3f ts;

    const Material *material;
    const TransportMode mode;

public:
    SeparableBSSRDF(
        const SurfaceInteraction &po,
        Float eta,
        const Material *material,
        TransportMode mode
    ) : BSSRDF(po, eta),
        ns(po.shading.n),
        ss(Normalize(po.shading.dpdu)),
        ts(Cross(ns, ss)),
        material(material),
        mode(mode)
    {}

    // Computes an approximation to the ratio of differential exitant radiance at po in
    // the direction SurfaceInteraction:wo to the incident differential flux at pi from
    // direction wi.
    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) {
        // TODO: this explanation doesn't make sense; I would have expected the indices
        // of refraction to be switched, so that the ray is exiting from within the surface
        // and into the outside.
        //
        // FrDielectric evaluates the Fresnel reflectance equation between 2 dielectric
        // media; the 1st argument is cosThetaI, the cosine of the angle of incidence
        // relative to the normal. 1 is the index of refraction of the medium the ray
        // is traveling in (1 corresponds to a vaccum) and eta is the index of refraction
        // of the medium within the surface that this BSSRDF represents.
        //
        // Ft is then the fraction that is transmitted into direction po.wo after exiting
        // the material.
        Float Ft = 1 - FrDielectric(Dot(po.wo, po.shading.n), 1, eta);

        // This is the separable approximation to the BSSRDF.
        return Ft * Sp(pi) * Sw(wi);
    }

    // Accounts for the influence of the boundary on the directional distribution of light
    // entering the surface from direction wi.
    Spectrum Sw(const Vector3f &wi) const {
        // c is a normalization that causes Sw(wi)*cos(Theta) to integrate to 1 over the
        // cosine-weighted hemisphere.
        Float c = 1 - 2 * FresnelMoment1(1 / eta);

        return (1 - FrDielectric(CosTheta(wi), 1, eta)) / (c * Pi);
    }

    // Gives the spatial distribution characterizing how far light travels after entering the
    // surface. PBRT doesn't give the actual expression of Sp and instead gives and implements
    // an approximation. This approximation is based on 2 assumptions:
    // 
	// - The surface is locally planar.
    //
    // - It isn't the location of the 2 points but the distance between them that affects the
    //   value of S.
    //
    // An auxiliary function Sr that is given the distance between po and pi gives the approximation
    // and is implemented by subclasses.
    Spectrum Sp(const SurfaceInteraction &pi) const {
        return Sr(Distance(po.p, pi.p));
    }

    virtual Spectrum Sr(Float d) const = 0;
};

// A tabulated BSSRDF is 5-dimensional (index of refraction eta, scattering anisotropy g, albedo rho,
// extinction coefficient sigma_t, and radius), but fixing sigma_t=1 and using optical radius 
// r_optical = r * sigma_t reduces the dimension. A Jacobian determinant factor is then needed to
// account for the radius-to-optical-radius change of variable.
struct BSSRDFTable {
    // The table's dimension is nRhoSamples * mOpticalRadiusSamples.

    // Number of recorded albedo samples.
    const int nRhoSamples;

    // Number of recorded optical radii samples. Optical radii r_optical = r * sigma_t is a unitless
    // transformation of the radius r = |po - pi|.
    const int mOpticalRadiusSamples;

    std::unique_ptr<Float[]> rhoSamples;
    std::unique_ptr<Float[]> opticalRadiusSamples;

    // Effective albedo, which is different from the single-scattering-event albedo stored in
    // TabulatedBSSRDF::rho. Effective albedo integrates the radial profile Sr over all radii r and
    // all polar directions, i.e. over an infinitely large disk? A nonlinear and strictly monotonically
    // increasing function of the single scattering albedo TabulatedBSSRDF::rho
    std::unique_ptr<Float[]> effectiveRho;

    // Stores a sample BSSRDF sample value for each of the
    // nRhoSamples * mOpticalRadiusSamples pairs (rho, r_optical).
    std::unique_ptr<Float[]> profile;
    std::unique_ptr<Float[]> profileCDF;

    BSSRDFTable(
        int nRhoSamples, 
        int mOpticalRadiusSamples
    ) : nRhoSamples(nRhoSamples),
        mOpticalRadiusSamples(mOpticalRadiusSamples),
        rhoSamples(new Float[nRhoSamples]),
        opticalRadiusSamples(new Float[mOpticalRadiusSamples]),
        profile(new Float[nRadiusSamples * nRhoSamples]),
        rhoEff(new Float[nRhoSamples]),
        profileCDF(new Float[nRadiusSamples * nRhoSamples])
    {}

    inline Float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * mOpticalRadiusSamples + radiusIndex];
    }
};

class TabulatedBSSRDF : public SeparableBSSRDF {
private:
    // A 2-dimensional, unitless table.
    const BSSRDFTable &table;

    // Extinction coefficient that gives the rate of scattering or absorption per unit distance.
    Spectrum sigma_t;

    // Albedo that gives the reduction in energy after a single scattering event.
    Spectrum rho;

public:
    TabulatedBSSRDF(
        const SurfaceInteraction &po,
        const Material *material,
        TransportMode mode,
        // Index of refraction.
        Float eta,
        const Spectrum &sigma_a,
        const Spectrum &sigma_s,
        const BSSRDFTable &table
    ) : SeparableBSSRDF(po, eta, material, mode), table(table) {
        sigma_t = sigma_a + sigma_s;
        for (int c = 0; c < Spectrum::nSamples; ++c) {
            rho[c] = sigma_t[c] != 0 ? (sigma_s[c] / sigma_t[c]) : 0;
        }
    }

    // Computes the radial profile of BSSRDF using spline-based interpolation of the tabulated
    // samples of r and rho.
    Spectrum Sr(Float r) const;
};

#endif // CPBRT_CORE_BSSRDF_H