#include "gaussian.h"
#include "core/paramset.h"

Float GaussianFilter::Evaluate(const Point2f &p) const {
    return Gaussian(p.x, gaussianAtRadiusX) * Gaussian(p.y, gaussianAtRadiusY);
}