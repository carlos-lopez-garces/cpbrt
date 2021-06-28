#include "mitchell.h"
#include "core/paramset.h"

Float MitchellFilter::Evaluate(const Point2f &p) const {
    return Mitchell1D(p.x * invRadius.x) * Mitchell1D(p.y * invRadius.y);
}

MitchellFilter *CreateMitchellFilter(const ParamSet &ps) {
    Float xw = ps.FindOneFloat("xwidth", 2.f);
    Float yw = ps.FindOneFloat("ywidth", 2.f);
    Float B = ps.FindOneFloat("B", 1.f / 3.f);
    Float C = ps.FindOneFloat("C", 1.f / 3.f);
    return new MitchellFilter(Vector2f(xw, yw), B, C);
}