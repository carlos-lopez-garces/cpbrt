#include "sampler.h"
#include "scene.h"
#include "spectrum.h"
#include "transform.h"

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
protected:

public:
    // Flags that characterize the light source.
    const int flags;
    
    const Transform LightToWorld;
    const Transform WorldToLight;

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

    // Computes incident radiance at the point of Interaction, as well as the direction
    // of incidence wi.
    virtual Spectrum Sample_Li(
        const Interaction &it,
        // This is not the point of incidence, but a sample point on the surface of the
        // light source. For "sampled" light sources that illuminate a given point from
        // different directions.
        const Point2f &u,
        Vector3f *wi,
        // Probability density function from sampled light sources.
        Float *pdf,
        // Occlusion information between the light source and the point of incidence.
        VisibilityTester *vis
    ) const = 0;

    // Computes total emitted power.
    virtual Spectrum Power() const = 0;

    // To be called during Scene construction and before rendering begins.
    virtual void Preprocess(const Scene &scene) {}
};

class VisibilityTester {
private:
    // Departure and destination points of a shadow ray.
    Interaction p0;
    Interaction p1;

public:
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