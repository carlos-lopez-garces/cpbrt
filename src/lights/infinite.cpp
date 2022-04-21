#include "lights/infinite.h"
#include "core/geometry.h"
#include "core/imageio.h"
#include "core/medium.h"

InfiniteAreaLight::InfiniteAreaLight(
    const Transform &lightToWorld,
    const Spectrum &power,
    int nSamples,
    const std::string &textureFilepath
) : Light((int) LightFlags::Infinite, lightToWorld, MediumInterface(), nSamples) {

    Point2i resolution;
    std::unique_ptr<RGBSpectrum[]> texels(nullptr);

    if (textureFilepath != "") {
        texels = ReadImage(textureFilepath, &resolution);
        if (texels) {
            for (int i = 0; i < resolution.x * resolution.y; ++i) {
                texels[i] *= power.ToRGBSpectrum();
            }
        }
    }

    // std::unique_ptr implements std::unique_ptr::operator bool, which returns
    // false when the unique_ptr is empty.
    if (!texels) {
        // No texture backs the light. Let the light be backed by a single-texel
        // texture with power for its value.
        resolution.x = resolution.y = 1;
        texels = std::unique_ptr<RGBSpectrum[]>(new RGBSpectrum[1]);
        texels[0] = power.ToRGBSpectrum();
    }

    // Create a mip map for the texture.
    LMap.reset(new MIPMap<RGBSpectrum>(resolution, texels.get()));

    // Initialize sampling PDFs.

    // Filter and scale the environment map to obtain scalar-valued image storing
    // the (filtered) luminance of texels that surround a sample point.
    int width = resolution.x;
    int height = resolution.y;
    // Filter size. Number/fraction of texels to filter per sample point.
    Float filterSize = (Float) 1 / std::max(width, height);
    std::unique_ptr<Float[]> scalarValuedImg(new Float[width * height]);
    for (int v = 0; v < height; ++v) {
        // Process vth row.
        Float vp = (Float) v / (Float) height;

        // Theta is a normalized fraction of Pi: kTheta, where k ~in (0.0,1.0).
        // An azimuth angle.
        Float sinTheta = std::sin(Pi * Float(v + 0.5f) / Float(height));

        for (int u = 0; u < width; ++u) {
            // Process uth column.
            Float up = (Float) u / (Float) width;
            scalarValuedImg[u + v * width] = LMap->Lookup(Point2f(up, vp), filterSize).y();
            // Multiplying by sinTheta is supposed to correct the distortion caused
            // by mapping the rectangular environment map to the interior surface of
            // the unit sphere.
            scalarValuedImg[u + v * width] *= sinTheta;
        }
    }

    // Compute sampling distributions for rows and columns of the image.
    distribution.reset(new Distribution2D(scalarValuedImg.get(), width, height));
}

Spectrum InfiniteAreaLight::Power() const {
    return Pi * worldRadius * worldRadius * Spectrum(LMap->Lookup(Point2f(.5f, .5f), .5f), SpectrumType::Illuminant);
}

Spectrum InfiniteAreaLight::Le(const RayDifferential &rd) const {
    Vector3f w = Normalize(WorldToLight(rd.d));
    // Texture sample point. Convert ray direction (a rectangular 3D coordinate)
    // to spherical coordinate.
    Point2f st(SphericalPhi(w) * Inv2Pi, SphericalTheta(w) * InvPi);
    return Spectrum(LMap->Lookup(st), SpectrumType::Illuminant);
}

Spectrum InfiniteAreaLight::Sample_Li(
    // A point on a surface possibly lit by this light.
    const Interaction &it,
    // Uniformy distributed 2D sample.
    const Point2f &u,
    // Sampled incident direction.
    Vector3f *wi,
    // Probability of sampling the returned wi.
    Float *pdf,
    VisibilityTester *vis
) {
    // Find (u,v) sample coordinates in environment map.
    Float mapPdf;
    Point2f uv = distribution->SampleContinuous(u, &mapPdf);
    if (mapPdf == 0) {
        return Spectrum(0.f);
    }

    // Turn infinite light sample point into spherical coordinates direction:
    // (theta, phi) = (vPi, 2uPi). Then turn spherical coordinate into rectangular
    // direction w = (x,y,z).
    Float theta = uv[1] * Pi;
    Float phi = 2* uv[0] * Pi;
    Float cosTheta = std::cos(theta);
    Float sinTheta = std::sin(theta);
    Float sinPhi = std::sin(phi);
    Float cosPhi = std::cos(phi);
    *wi = LightToWorld(Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));

    // Compute PDF for sampled direction. The sampled (u,v) underwent 2 mappings:
    // g: (u,v) -> (theta, phi)
    // h: (theta, phi) -> (x,y,z) (latitude-longitude mapping)
    //
    // As is the case with any change of variables, we need to multiply by a scaling
    // factor (by the determinant of the Jacobian matrix, if multivariable. The
    // determinant of the Jacobian of g is |Jg|=2Pi^2. The determinant of the Jacobian
    // of h is |Jh|=sin(theta). Thus, the scaling factor of the mapping gh: (u,v) -> (x,y,z)
    // is |Jgh|=2sin(theta)Pi^2.
    if (sinTheta == 0) {
        *pdf = 0;
    } else {
        *pdf = mapPdf / (2 * Pi * Pi * sinTheta);
    }

    // Return radiance along sampled direction. We mapped the environment map to the interior
    // of the unit sphere. However, the light sample must come from infinitely far away and
    // any or no object may be in the way between the Interaction point and the light sample,
    // so we construct a point beyond the bounds of the scene and along the sampled direction
    // of incidence to test for visibility.
    *vis = VisibilityTester(it, Interaction(it.p + *wi * (2 * worldRadius), it.time, mediumInterface));
    return Spectrum(LMap->Lookup(uv), SpectrumType::Illuminant);
}

Float InfiniteAreaLight::Pdf_Li(const Interaction &it, const Vector3f &w) const {
    Vector3f wi = WorldToLight(w);
    Float theta = SphericalTheta(wi);
    Float phi = SphericalPhi(wi);
    Float sinTheta = std::sin(theta);
    if (sinTheta == 0) {
        return 0;
    }

    // Apply the inverse of the mapping g: (theta, phi) -> (vPi, 2uPi) to recover the (u,v) coordinate
    // sampled by Sample_Li from the environment map and that corresponds to the input, world-space
    // direction w:
    // g^-1: (u,v) -> (phi/2Pi, theta/Pi)
    //
    // Then, scale p(u,v) using the Jacobian determinants of g and the latitude-longitude mapping
    // h: (theta, phi) -> (x,y,z).  
    return distribution->Pdf(Point2f(phi * Inv2Pi, theta * InvPi)) / (2 * Pi * Pi * sinTheta);
}