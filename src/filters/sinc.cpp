#include "sinc.h"
#include "core/paramset.h"

Float LanczosSincFilter::Evaluate(const Point2f &p) const {
    return WindowedSinc(p.x, radius.x) * WindowedSinc(p.y, radius.y);
}

LanczosSincFilter *CreateSincFilter(const ParamSet &ps) {
    Float xw = ps.FindOneFloat("xwidth", 4.);
    Float yw = ps.FindOneFloat("ywidth", 4.);
    Float tau = ps.FindOneFloat("tau", 3.f);
    return new LanczosSincFilter(Vector2f(xw, yw), tau);
}