#include "diffuse.h"

Spectrum DiffuseAreaLight::Power() const {
    return (twoSided ? 2 : 1) * emittedL * area * Pi;
}

Spectrum DiffuseAreaLight::Sample_Li(
    const Interaction &it,
    const Point2f &u,
    Vector3f *wi,
    Float *pdf,
    VisibilityTester *vis
) const {
    // Sample a point on the surface of the area light.
    Interaction sampleSurfacePoint = shape->Sample(it, u, pdf);
    
    sampleSurfacePoint.mediumInterface = mediumInterface;
    
    // The sampled direction of incidence is the normalized vector that goes from
    // the sampled point on the surface of the light source to the point of intersection
    // on the illuminated object's surface.
    *wi = Normalize(sampleSurfacePoint.p - it.p);

    // Corresponding sampling probability.
    *pdf = shape->Pdf(it, *wi);

    // Is the light source occluded in the sampled direction?
    *vis = VisibilityTester(it, sampleSurfacePoint);

    // Evaluate emitted radiance.
    return L(sampleSurfacePoint, -*wi);
}

Float DiffuseAreaLight::Pdf_Li(const Interaction &it, const Vector3f &wi) const {
    return shape->Pdf(it, wi);
}

std::shared_ptr<AreaLight> CreateDiffuseAreaLight(
    const Transform &LightToWorld, 
    const Medium *medium,
    const ParamSet &params,
    const std::shared_ptr<Shape> &shape
) {
    Spectrum L = params.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum scale = params.FindOneSpectrum("scale", Spectrum(1.0));
    int nSamples = params.FindOneInt("samples", params.FindOneInt("nsamples", 1));
    bool twoSided = params.FindOneBool("twosided", false);
    return std::make_shared<DiffuseAreaLight>(LightToWorld, medium, L * scale, nSamples, shape, twoSided);
}