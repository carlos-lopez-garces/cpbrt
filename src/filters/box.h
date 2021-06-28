#ifndef CPBRT_FILTERS_BOX_H
#define CPBRT_FILTERS_BOX_H

#include "core/filter.h"

class BoxFilter : public Filter {
public:
    BoxFilter(const Vector2f &radius) : Filter(radius) {}

    // The box filter maps any sample point to a constant value. Sample point coordinates
    // are expressed with respect to the filter's center and normalized to [-radius, radius],
    // the filter's extent.
    Float Evaluate(const Point2f &p) const;
};

BoxFilter *CreateBoxFilter(const ParamSet &ps);

#endif // CPBRT_FILTERS_BOX_H