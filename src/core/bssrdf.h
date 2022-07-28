#ifndef CPBRT_CORE_BSSRDF_H
#define CPBRT_CORE_BSSRDF_H

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

#endif // CPBRT_CORE_BSSRDF_H