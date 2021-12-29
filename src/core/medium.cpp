#include "medium.h"
#include "geometry.h"

Float HenyeyGreensteinPhaseFunction::p(const Vector3f &wo, const Vector3f &wi) const {
    return PhaseHG(Dot(wo, wi), g);
}