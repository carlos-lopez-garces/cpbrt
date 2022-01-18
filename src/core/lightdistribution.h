#ifndef CPBRT_CORE_LIGHTDISTRIBUTION_H
#define CPBRT_CORE_LIGHTDISTRIBUTION_H

#include "cpbrt.h"
#include "geometry.h"
#include "sampling.h"
#include "scene.h"

class LightDistribution {
public:
    virtual ~LightDistribution() {}

    virtual const Distribution1D *Lookup(const Point3f &p) const = 0;
};

class UniformLightDistribution : public LightDistribution {
public:
    UniformLightDistribution(const Scene &scene);

    const Distribution1D *Lookup(const Point3f &p) const;

private:
    std::unique_ptr<Distribution1D> distrib;
};

#endif // CPBRT_CORE_LIGHTDISTRIBUTION_H