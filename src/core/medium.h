#ifndef PBRT_CORE_MEDIUM_H
#define PBRT_CORE_MEDIUM_H

#include "cpbrt.h"
#include "geometry.h"
#include "spectrum.h"

struct MediumInterface {
    // TODO: implement.
    const Medium *inside;
    const Medium *outside;

    MediumInterface() : inside(nullptr), outside(nullptr) {}
    
    MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}
    
    MediumInterface(const Medium *inside, const Medium *outside)
        : inside(inside), outside(outside) {}
};

#endif // PBRT_CORE_MEDIUM_H