#ifndef CPBRT_CORE_LIGHT_H
#define CPBRT_CORE_LIGHT_H

#include "cpbrt.h"
#include "memory.h"
#include "interaction.h"

enum class LightFlags : int {
    // Single position delta distribution.
    DeltaPosition = 1,

    // Single direction delta distribution.
    DeltaDirection = 2,

    // Area light.
    Area = 4,

    // Infinite area light.
    Infinite = 8
};

inline bool IsDeltaLight(int flags) {
    return flags & (int) LightFlags::DeltaPosition
           || flags & (int) LightFlags::DeltaDirection;
}

class Light {
public:
    // Flags that characterize the light source.
    const int flags;
    
    const Transform LightToWorld;
    const Transform WorldToLight;

    // Number of (point) samples that an Integrator should take from this light source.
    const int nSamples;

    // Participating medium on the inside and the outside of the light source.
    // nullptr for vacuum.
    const MediumInterface mediumInterface;

    Light(
        int flags,
        const Transform &LightToWorld,
        const MediumInterface &mediumInterface,
        int nSamples = 1
    ) : flags(flags),
        nSamples(std::max(1, nSamples)),
        mediumInterface(mediumInterface),
        LightToWorld(LightToWorld),
        WorldToLight(Inverse(LightToWorld)) {
        
        // TODO: Warn if light has transformation with non-uniform scale.
    }

    // Samples the direction of incidence wi at the point of Interaction and computes
    // the corresponding incident radiance.
    virtual Spectrum Sample_Li(
        const Interaction &it,
        // This is not the point of incidence, but a 2D random sample value used for
        // sampling this light source. For light sources that illuminate a given point
        // from different directions (like those that have a surface).
        const Point2f &u,
        // Sampled direction.
        Vector3f *wi,
        // Probability density function of distribution of directions that correspond
        // to this light source.
        Float *pdf,
        // Occlusion information between the light source and the point of incidence.
        VisibilityTester *vis
    ) const = 0;

    // Evaluates the PDF of the distribution of directions of this light source for
    // the input incident direction.
    virtual Float Pdf_Li(const Interaction &it, const Vector3f &wi) const = 0;

    // Computes total emitted power.
    virtual Spectrum Power() const = 0;

    // Computes emitted radiance by this light source in the direction of the given ray.
    virtual Spectrum Le(const RayDifferential &rd) const;

    // To be called during Scene construction and before rendering begins.
    virtual void Preprocess(const Scene &scene) {}
};

class VisibilityTester {
private:
    // Departure and destination points of a shadow ray.
    Interaction p0;
    Interaction p1;

public:
    VisibilityTester() {}

    VisibilityTester(const Interaction &p0, const Interaction &p1) : p0(p0), p1(p1) {}

    const Interaction &P0() const {
        return p0;
    }

    const Interaction &P1() const {
        return p1;
    }

    // Tells whether any occluder is intersected by a shadow ray between p0 and p1,
    // ignoring participating media that may have an effect on the ray's radiance.
    bool Unoccluded(const Scene &scene) const;

    // TODO: implement when implementing volume scattering.
    Spectrum Tr(const Scene &scene, Sampler &sampler) const;
};

#endif // CPBRT_CORE_LIGHT_H