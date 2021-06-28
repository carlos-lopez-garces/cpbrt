#ifndef CPBRT_FILTERS_SINC_H
#define CPBRT_FILTERS_SINC_H

#include <cmath>

#include "core/filter.h"

class LanczosSincFilter : public Filter {
private:
    // Controls the number of cycles of the windowed sinc function.
    const Float tau;

public:
    LanczosSincFilter(const Vector2f &radius, Float tau)
        : Filter(radius), tau(tau)
    {}

    Float Evaluate(const Point2f &p) const;

    // The sinc filter has infinite extent and doesn't fall off to 0. To reconstruct
    // correctly, it needs an infinite number of samples, which isn't possible in rendering.
    Float Sinc(Float sampleCoordinate) const {
        Float x = sampleCoordinate;

        x = std::abs(x);
        if (x < 1e-5) {
            return 1.f;
        }
        return std::sin(Pi * x) / (Pi * x);
    }

    // The Lanczos windowed sinc function limits the extent of the sinc function to a
    // number of cycles controlled by the tau parameter. It also makes it fall off to 0
    // at the radius.
    Float WindowedSinc(Float sampleCoordinate, Float radius) const {
        Float x = sampleCoordinate;

        x = std::abs(x);
        if (x > radius) {
            return 0.f;
        }
        Float lanczos = Sinc(x / tau);
        return Sinc(x) * lanczos;
    }
};

LanczosSincFilter *CreateSincFilter(const ParamSet &ps);

#endif // PBRT_FILTERS_SINC_H