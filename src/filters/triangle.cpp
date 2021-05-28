#include <algorithm>
#include <cmath>

#include "triangle.h"
#include "core/paramset.h"

Float TriangleFilter::Evaluate(const Point2f &p) const {
    return std::max((Float) 0, radius.x - std::abs(p.x)) 
           * std::max((Float) 0, radius.y - std::abs(p.y));
}