#ifndef CPBRT_FILTERS_GAUSSIAN_H
#define CPBRT_FILTERS_GAUSSIAN_H

#include <cmath>

#include "core/filter.h"

class GaussianFilter : public Filter {
private:
    // Alpha determines the rate of decay of the curve.
    const Float alpha;
    const Float gaussianAtRadiusX;
    const Float gaussianAtRadiusY;

public:
    GaussianFilter(const Vector2f &radius, Float alpha)
        : Filter(radius),
          alpha(alpha),
          gaussianAtRadiusX(std::exp(-alpha * radius.x * radius.x)),
          gaussianAtRadiusY(std::exp(-alpha * radius.y * radius.y))
    {}

    // The Gaussian filter maps sample points to values of a bell curve that decrease
    // nonlinearly as the coordinates of the sample point tend to the filter's radius,
    // at which point the function's value is 0.
    Float Evaluate(const Point2f &p) const;

    Float Gaussian(Float sampleCoordinate, Float gaussianAtRadius) const {
        // Subtracting the Gaussian value at the radius from the Gaussian value at the 
        // sample point coordinate ensures that it tends and reaches 0 as the coordinate
        // tends to and reaches the radius.
        return std::max(
            (Float) 0,
            Float(std::exp(-alpha * sampleCoordinate * sampleCoordinate) - gaussianAtRadius)
        );
    }
};

#endif // CPBRT_FILTERS_GAUSSIAN_H