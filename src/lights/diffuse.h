#ifndef CPBRT_LIGHTS_DIFFUSE_H
#define CPBRT_LIGHTS_DIFFUSE_H

#include "core/cpbrt.h"
#include "core/light.h"
#include "core/paramset.h"
#include "core/primitive.h"

class DiffuseAreaLight : public AreaLight {
protected:
    // Emitted radiance (constant across the surface).
    const Spectrum emittedL;

    std::shared_ptr<Shape> shape;
    const Float area;

    // Whether the light source emits light from both its front and back faces.
    const bool twoSided;

public:
    DiffuseAreaLight(
        const Transform &LightToWorld,
        const MediumInterface &mediumInterface,
        const Spectrum &emittedL,
        int nSamples,
        const std::shared_ptr<Shape> &shape,
        bool twoSided = false
    ) : AreaLight(LightToWorld, mediumInterface, nSamples),
        emittedL(emittedL),
        shape(shape),
        area(shape->Area()),
        twoSided(twoSided)
    {}

    Spectrum L(const Interaction &intr, const Vector3f &w) const {
        // Only emit light from the face of the Shape where the normal is outward-facing.
        // If the ray intersects the back face, no radiance is returned.
        return (twoSided || Dot(intr.n, w) > 0.0f) ? emittedL : Spectrum(0.0f);
    }

    Spectrum Power() const;

    // Because the Interaction point 'it' on the surface of the intersected object may be illuminated by
    // this area light from multiple directions, Sample_Li samples an incident direction wi uniformly at
    // random from the directional distribution of the area light (an incident direction that can then be
    // used to sample the intersected object's BxDF). It also computes the radiance emitted by the light
    // source in that direction.
    virtual Spectrum Sample_Li(
        // Intersection point on the surface of an object that may be illuminated by this light source
        // (unless it is occluded).
        const Interaction &it,
        // A 2D uniform random sample used for sampling the directional distribution of the light source.
        const Point2f &u,
        // The returned sampled direction.
        Vector3f *wi,
        // Probability with which the direction was sampled, obtained from the density function of the
        // directional distribution of the light source.
        Float *pdf,
        // Occlusion information between the light source and the point of incidence.
        VisibilityTester *vis
    ) const;

    Float Pdf_Li(const Interaction &it, const Vector3f &wi) const;
};

std::shared_ptr<AreaLight> CreateDiffuseAreaLight(
    const Transform &LightToWorld, 
    const Medium *medium,
    const ParamSet &params,
    const std::shared_ptr<Shape> &shape
);

#endif // CPBRT_LIGHTS_DIFFUSE_H