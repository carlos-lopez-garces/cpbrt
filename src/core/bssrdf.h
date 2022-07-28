#ifndef CPBRT_CORE_BSSRDF_H
#define CPBRT_CORE_BSSRDF_H

#include "geometry.h"
#include "interaction.h"

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

// A separable BSSRDF approximates a BSSRDF with a product of 3 functions:
//
// S(po, wo, pi, wi) ~= (1 - Fr(cos(Theta_o)))Sp(po,pi)Sw(wi)
//
// where Fr is Fresnel transmittance and 1 - Fr(cos(Theta_o)) is the fraction of incident
// radiance that gets transmitted into direction wo when exiting the surface; Sw(wi) is
// the fraction of incident radiance coming from wi and that makes it from pi to po as a 
// result of the influence of the characteristics of the surface boundary at pi; and Sp(po,pi)
// is a spatial distribution that tells how far light travels within the surface.
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
};

#endif // CPBRT_CORE_BSSRDF_H