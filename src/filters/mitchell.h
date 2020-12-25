#include <cmath>

#include "core/filter.h"

class MitchellFilter : public Filter {
private:
    // B and C are chosen to satisfy B + 2C = 1 and so that the Mitchell-Netravali
    // polynomial is continuous on [-2, 2]. Bad values of B and C may cause jump
    // discontinuities because the function is defined piecewise over the [0, 1] and
    // [1, 2] subdomains (and [-2, -1] and [-1, 0] by reflection about the y-axis).
    const Float B;
    const Float C;

public:
    MitchellFilter(const Vector2f &radius, Float B, Float C) 
        : Filter(radius), B(B), C(C)
    {}

    Float Evaluate(const Point2f &p) const;

    // The Mitchell-Netravali 1D function. The function is separable, so computing it
    // in 2D amounts to multiplying the function in 1D of the 2 dimensions.
    //
    // The function is a piecewise defined cubic polynomial with subdomains [-2, -1],
    // [-1, 0], [0, 1], and [1, 2].
    //
    // The function is negative by design at some trailing subinterval. This helps give
    // sharpness to edges.
    Float Mitchell1D(Float sampleCoordinate) const {
        Float x = sampleCoordinate;

        x = std::abs(2 * x);
        if (x > 1) {
            return ((-B - 6*C)*x*x*x + (6*B + 30*C)*x*x + (-12*B - 48*C)*x + (8*B + 24*C)) * (1.f/6.f);
        } else {
            return ((12 - 9*B - 6*C)*x*x*x + (-18 + 12*B + 6*C)*x*x + (6 - 2*B)) * (1.f/6.f);
        }
    }
};