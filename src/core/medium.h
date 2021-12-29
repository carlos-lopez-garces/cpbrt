#ifndef PBRT_CORE_MEDIUM_H
#define PBRT_CORE_MEDIUM_H

#include "cpbrt.h"
#include "geometry.h"
#include "spectrum.h"

struct MediumInterface {
    // TODO: implement.
    const Medium *inside;
    const Medium *outside;

    MediumInterface() : inside(nullptr), outside(nullptr) {}
    
    MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}
    
    MediumInterface(const Medium *inside, const Medium *outside)
        : inside(inside), outside(outside) {}
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