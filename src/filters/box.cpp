#include "box.h"
#include "core/paramset.h"

Float BoxFilter::Evaluate(const Point2f &p) const {
    return 1.0;
}

BoxFilter *CreateBoxFilter(const ParamSet &ps) {
    Float xw = ps.FindOneFloat("xwidth", 0.5f);
    Float yw = ps.FindOneFloat("ywidth", 0.5f);
    return new BoxFilter(Vector2f(xw, yw));
}