#include "light.h"

bool VisibilityTester::Unoccluded(const Scene &scene) const {
    // TODO: implement Interaction::SpawnRayTo.
    return !scene.IntersectP(p0.SpawnRayTo(p1));
}