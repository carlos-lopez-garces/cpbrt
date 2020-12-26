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