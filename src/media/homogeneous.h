#ifndef CPBRT_MEDIA_HOMOGENEOUS_H
#define CPBRT_MEDIA_HOMOGENEOUS_H

#include "core/medium.h"
#include "core/spectrum.h"

class HomogeneousMedium : public Medium {
private:
    // Absorption cross section. The probability that radiance is absorbed per unit
    // distance. Absorption is a source of attenuation.
    const Spectrum sigma_a;

    // Scattering coefficient. The probability of occurrence of an out-scattering event
    // per unit distance, that is, the probability that radiance exits a differential
    // region of the beam due to deflection by collision with particles. Out-scattering
    // is a source of attenuation.
    const Spectrum sigma_s;

    // Attenuation or extinction coefficient. The combined effect of absorption and
    // outscattering: sigma_t(p,w) = sigma_a(p,w) + sigma_s(p,w). 
    const Spectrum sigma_t;

    // Henyey-Greenstein phase function's asymmetry parameter. Negative g values describe
    // phase functions that primarily scatter light back in the incident direction, and
    // positive g values, forward in the direction that light is travelling. Isotropic
    // phase functions, which describe equal scattering in all directions, have g = 0.
    const Float g;

public:
    // These coefficients are constant across the extent of an homogeneous medium.
    HomogeneousMedium(
        const Spectrum &sigma_a,
        const Spectrum &sigma_s,
        Float g
    ) : sigma_a(sigma_a), sigma_s(sigma_s), sigma_t(sigma_a + sigma_s), g(g) {}

    // Computes the beam transmittance, the fraction of radiance that is transmitted 
    // through the medium between 2 points, in this case, between the origin of the
    // ray and the point that corresponds to tMax.
    //
    // The beam transmittance equation is Tr(p->p') = e^(-tau(p->p')), where
    // tau(p->p') = sigma_t*d and d = ||p - p'|| according to Beer's law, because
    // sigma_a and sigma_s are constant across the medium, so they don't need to be
    // integrated with respect to p along p - p'.
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;

    Spectrum Sample(
        const Ray &ray,
        Sampler &sampler,
        MemoryArena &arena,
        MediumInteraction *mi
    ) const;
};

#endif // CPBRT_MEDIA_HOMOGENEOUS_H