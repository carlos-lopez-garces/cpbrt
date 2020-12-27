#include <cmath>

#include "film.h"

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

        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) {
            mergePixel.xyz[i] += xyz[i];
        }

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