#ifndef PBRT_CORE_MEDIUM_H
#define PBRT_CORE_MEDIUM_H

#include "cpbrt.h"
#include "geometry.h"
#include "spectrum.h"

// MediumInterfaces are associated with GeometricPrimitives to define the
// scattering media that exists inside and outside of them. Null pointers
// indicate a vacuum.
struct MediumInterface {
    const Medium *inside;
    const Medium *outside;

    MediumInterface() : inside(nullptr), outside(nullptr) {}
    
    // Same medium on both sides.
    MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}
    
    MediumInterface(const Medium *inside, const Medium *outside)
        : inside(inside), outside(outside) {}
    
    // Checks whether the interface marks a transition between 2 different media.
    bool IsMediumTransition() const {
        return inside != outside;
    }

    // Samples the medium's scattering interaction along the ray.
    virtual Spectrum Sample(
        const Ray &ray,
        Sampler &sampler,
        MemoryArena &arena,
        MediumInteraction *mi
    ) const = 0;
};

class Medium {
public:
    virtual ~Medium() {}

    // Computes the beam transmittance, the fraction of radiance that is transmitted 
    // through the medium between 2 points, in this case, between the origin of the
    // ray and the point that corresponds to tMax. Implementations must account for
    // the effects of absorption and out-scattering.
    //
    // The beam transmittance equation is Tr(p->p') = e^(-tau(p->p')), where
    // tau(p->p') = \int{0}{d}{sigma_t(p + tw, w)}dt and d = ||p - p'||.
    //
    // The input ray is assumed to be fully contained by the medium and that there
    // aren't occluders between its origin and the point that corresponds to tMax.
    virtual Spectrum Tr(const Ray &ray, Sampler &sampler) const = 0;
    
    virtual Spectrum Sample(
        const Ray &ray,
        Sampler &sampler,
        MemoryArena &arena,
        MediumInteraction *mi
    ) const = 0;
};

class PhaseFunction {
public:
    // Evaluates the phase function for the given incident and outgoing directions.
    // Phase functions are reciprocal, so p(wo, wi) = p(wi, wo).
    virtual Float p(const Vector3f &wo, const Vector3f &wi) const = 0;
    // Sample the phase function, choosing wi probabilistically.
    virtual Float Sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const = 0;
};

// Henyey-Greenstein phase function, designed to be easy to fit to measured scattering
// data. The distribution of scattering light is controlled by the asymmetry parameter g.
inline Float PhaseHG(Float cosTheta, Float g) {
    Float denominator = 1 + g*g + 2*g*cosTheta;
    return Inv4Pi * (1 - g*g) / (denominator * std::sqrt(denominator));
}

class HenyeyGreensteinPhaseFunction : public PhaseFunction {
private:
    // Asymmetry parameter. Controls the directional distribution of scattering. If negative,
    // light scatters back primarily; if positive, light scatters forward in the incident
    // direction.
    const Float g;

public:
    HenyeyGreensteinPhaseFunction(Float g) : g(g) {}

    // Evaluates the Henyey-Greenstein phase function for the given incident and outgoing
    // directions. Phase functions are reciprocal, so p(wo, wi) = p(wi, wo).
    virtual Float p(const Vector3f &wo, const Vector3f &wi) const = 0;

    // Sample the phase function, choosing wi probabilistically.
    virtual Float Sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const = 0;
};

#endif // PBRT_CORE_MEDIUM_H