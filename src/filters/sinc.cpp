#include "sinc.h"

Float LanczosSincFilter::Evaluate(const Point2f &p) const {
    return WindowedSinc(p.x, radius.x) * WindowedSinc(p.y, radius.y);
}