#ifndef CPBRT_CORE_FILTER_H
#define CPBRT_CORE_FILTER_H

#include "cpbrt.h"

class Filter {
public:
    const Vector2f radius;
    const Vector2f invRadius;

    // 2*radius is the filter's extent, i.e. the interval that it spans before becoming
    // or falling off to 0.
    Filter(const Vector2f &radius) 
        : radius(radius), invRadius(Vector2f(1 / radius.x, 1 / radius.y)) {}

    // Computes the weight for a sample point (sampleX, sampleY) expressed in
    // relation to the filter's origin (x, y), that is, (x - sampleX, y - sampleY).
    virtual Float Evaluate(const Point2f &p) const = 0;
};

#endif // CPBRT_CORE_FILTER_H