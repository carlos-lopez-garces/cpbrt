#ifndef CPBRT_CORE_LIGHTDISTRIBUTION_H
#define CPBRT_CORE_LIGHTDISTRIBUTION_H

#include "cpbrt.h"
#include "geometry.h"

class LightDistribution {
  public:
    virtual ~LightDistribution() {}

    virtual const Distribution1D *Lookup(const Point3f &p) const = 0;
};

#endif // CPBRT_CORE_LIGHTDISTRIBUTION_H