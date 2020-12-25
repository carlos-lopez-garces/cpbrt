#include "gaussian.h"

Float GaussianFilter::Evaluate(const Point2f &p) const {
    return Gaussian(p.x, gaussianAtRadiusX) * Gaussian(p.y, gaussianAtRadiusY);
}