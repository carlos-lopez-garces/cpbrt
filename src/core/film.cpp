#include <cmath>

#include "api.h"
#include "film.h"
#include "paramset.h"
#include "imageio.h"

Bounds2i Film::GetSampleBounds() const {
    Bounds2f floatBounds(
        // Offset the croppedPixelBounds outward to accomodate for the extent of the filters
        // centered at boundary pixels.
        Floor(Point2f(croppedPixelBounds.pMin) + Vector2f(0.5f, 0.5f) - filter->radius),
        Ceil(Point2f(croppedPixelBounds.pMax) - Vector2f(0.5f, 0.5f) + filter->radius)
    );

    return (Bounds2i) floatBounds;
}

Bounds2f Film::GetPhysicalExtent() const {
    Float aspect = (Float) fullResolution.y / (Float) fullResolution.x;
    // Pythagorean relation.
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;

    return Bounds2f(Point2f(-x/2, -y/2), Point2f(x/2, y/2));
}

// sampleBounds is the discrete pixel area that the film tile will cover. By handing out
// a std::unique_ptr to the FilmTile to the caller, the function is transferring ownership of
// the tile to the caller.
std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in sampleBounds contribute to.
    
    // Discrete-to-continuous pixel coordinate mapping is: c = d - 0.5.
    // Conversion to continuous coordinates is necessary because the filter's radius is a
    // floating-point number.
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f) sampleBounds;
    // Upper-left corner of discrete pixel area expanded by filter's radius. (The y-axis increases downwards.)
    Point2i p0 = (Point2i) Ceil(floatBounds.pMin - halfPixel - filter->radius);
    // Lower-right corner.
    Point2i p1 = (Point2i) Floor(floatBounds.pMax - halfPixel + filter->radius) + Point2i(1, 1);

    // Shave off the pixels that are outside the image.
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);

    return std::unique_ptr<FilmTile>(
        new FilmTile(tilePixelBounds, filter->radius, filterTable, filterTableWidth)
    );
}

// By receiving the FilmTile in a std::unique_ptr, the function is receiving ownership of the
// tile from the caller. Upon return, when the std::unique_ptr goes out of scope, the tile gets
// freed.
void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge pixel into Film::pixels.

        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);

        // The numerator in the pixel reconstruction/filtering equation is a Spectrum that can
        // give its value in XYZ color.
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) {
            mergePixel.xyz[i] += xyz[i];
        }

        // The denominator in the pixel reconstruction/filtering equation.
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

// ?
void Film::SetImage(const Spectrum *image) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        image[i].ToXYZ(p.xyz);
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}

// ?
void Film::AddSplat(const Point2f &p, const Spectrum &v) {
    if (!InsideExclusive((Point2i)p, croppedPixelBounds)) {
        return;
    }
    
    Float xyz[3];
    v.ToXYZ(xyz);
    
    Pixel &pixel = GetPixel((Point2i)p);
    for (int i = 0; i < 3; ++i) {
        pixel.splatXYZ[i].Add(xyz[i]);
    }
}

void Film::WriteImage(Float splatScale) {
    // Convert image to RGB and compute final pixel values.
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB. The pixel's XYZ color comes from the numerator
        // of the pixel reconstruction/filtering equation, which is a Spectrum built up from
        // filtered sample contributions.
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[3 * offset]);

        // Complete the computation of the pixel reconstruction/filtering equation by dividing by
        // the sum of the sample filter weights.
        Float filterWeightSum = pixel.filterWeightSum;
        if (filterWeightSum != 0) {
            Float reciprocalFilterWeight = (Float) 1 / filterWeightSum;
            // Clamp to 0 because reciprocalFilterWeight may be negative because some filters have
            // negative intervals (MitchellFilter, for example, has negative lobes intended to give
            // sharpness to edges).
            rgb[3 * offset  ] = std::max((Float) 0, rgb[3 * offset  ] * reciprocalFilterWeight);
            rgb[3 * offset+1] = std::max((Float) 0, rgb[3 * offset+1] * reciprocalFilterWeight);
            rgb[3 * offset+2] = std::max((Float) 0, rgb[3 * offset+2] * reciprocalFilterWeight);
        }

        // TODO: Explain splatting.
        // Add splat value at pixel.
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1], pixel.splatXYZ[2]};
        XYZToRGB(splatXYZ, splatRGB);
        rgb[3 * offset  ] += splatScale * splatRGB[0];
        rgb[3 * offset+1] += splatScale * splatRGB[1];
        rgb[3 * offset+2] += splatScale * splatRGB[2];

        // Scale pixel value by user-supplied scale.
        rgb[3 * offset  ] *= scale;
        rgb[3 * offset+1] *= scale;
        rgb[3 * offset+2] *= scale;

        ++offset;
    }

    // Write RGB image.
    ::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}

// TODO: shares a lot of code with Film::WriteImage.
std::unique_ptr<uint8_t[]> Film::GetPixels(Float splatScale) {
    // Convert image to RGB and compute final pixel values.
    std::unique_ptr<Float[]> rgb(new Float[4 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB. The pixel's XYZ color comes from the numerator
        // of the pixel reconstruction/filtering equation, which is a Spectrum built up from
        // filtered sample contributions.
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[4 * offset]);

        // Complete the computation of the pixel reconstruction/filtering equation by dividing by
        // the sum of the sample filter weights.
        Float filterWeightSum = pixel.filterWeightSum;
        if (filterWeightSum != 0) {
            Float reciprocalFilterWeight = (Float) 1 / filterWeightSum;
            // Clamp to 0 because reciprocalFilterWeight may be negative because some filters have
            // negative intervals (MitchellFilter, for example, has negative lobes intended to give
            // sharpness to edges).
            rgb[4 * offset  ] = std::max((Float) 0, rgb[4 * offset  ] * reciprocalFilterWeight);
            rgb[4 * offset+1] = std::max((Float) 0, rgb[4 * offset+1] * reciprocalFilterWeight);
            rgb[4 * offset+2] = std::max((Float) 0, rgb[4 * offset+2] * reciprocalFilterWeight);
            rgb[4 * offset+3] = 1;
        }

        // TODO: Explain splatting.
        // Add splat value at pixel.
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1], pixel.splatXYZ[2]};
        XYZToRGB(splatXYZ, splatRGB);
        rgb[4 * offset  ] += splatScale * splatRGB[0];
        rgb[4 * offset+1] += splatScale * splatRGB[1];
        rgb[4 * offset+2] += splatScale * splatRGB[2];

        // Scale pixel value by user-supplied scale.
        rgb[4 * offset  ] *= scale;
        rgb[4 * offset+1] *= scale;
        rgb[4 * offset+2] *= scale;

        ++offset;
    }
    
    // 8-bit formats; apply gamma
    Vector2i resolution = croppedPixelBounds.Diagonal();
    std::unique_ptr<uint8_t[]> rgb8(
        new uint8_t[4 * resolution.x * resolution.y]);
    uint8_t *dst = rgb8.get();
    for (int y = 0; y < resolution.y; ++y) {
        for (int x = 0; x < resolution.x; ++x) {
#define TO_BYTE(v) (uint8_t) Clamp(255.f * GammaCorrect(v) + 0.5f, 0.f, 255.f)
            dst[0] = TO_BYTE(rgb[4 * (y * resolution.x + x) + 0]);
            dst[1] = TO_BYTE(rgb[4 * (y * resolution.x + x) + 1]);
            dst[2] = TO_BYTE(rgb[4 * (y * resolution.x + x) + 2]);
            dst[3] = 1;
#undef TO_BYTE
            dst += 4;
        }
    }

    return rgb8;
}

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    std::string filename;
    if (CpbrtOptions.imageFile != "") {
        filename = CpbrtOptions.imageFile;
        std::string paramsFilename = params.FindOneString("filename", "");
        if (paramsFilename != "") {
            Warning(
                "Output filename supplied on command line, \"%s\" is overriding "
                "filename provided in scene description file, \"%s\".",
                CpbrtOptions.imageFile.c_str(), paramsFilename.c_str());
        }
    } else {
        filename = params.FindOneString("filename", "pbrt.exr");
    }

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (CpbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (CpbrtOptions.quickRender) yres = std::max(1, yres / 4);

    Bounds2f crop;
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr) {
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);
    } else {
        crop = Bounds2f(Point2f(Clamp(CpbrtOptions.cropWindow[0][0], 0, 1),
                                Clamp(CpbrtOptions.cropWindow[1][0], 0, 1)),
                        Point2f(Clamp(CpbrtOptions.cropWindow[0][1], 0, 1),
                                Clamp(CpbrtOptions.cropWindow[1][1], 0, 1)));
    }

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    return new Film(
        Point2i(xres, yres), crop, std::move(filter), diagonal, filename, scale
    );
}