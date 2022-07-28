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
};

#endif // CPBRT_CORE_BSSRDF_H