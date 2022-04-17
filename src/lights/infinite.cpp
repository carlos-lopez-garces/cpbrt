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
}

Spectrum InfiniteAreaLight::Power() const {
    return Pi * worldRadius * worldRadius * Spectrum(LMap->Lookup(Point2f(.5f, .5f), .5f), SpectrumType::Illuminant);
}