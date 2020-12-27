#include "core/filter.h"

class TriangleFilter : public Filter {
public:
    TriangleFilter(const Vector2f &radius) : Filter(radius) {}

    // The triangle filter maps sample points to values that decrease/decay/fall-off
    // linearly as the distance of the point to the filter's center increases.
    Float Evaluate(const Point2f &p) const;
};