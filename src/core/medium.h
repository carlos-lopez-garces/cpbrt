#ifndef PBRT_CORE_MEDIUM_H
#define PBRT_CORE_MEDIUM_H

#include "cpbrt.h"

struct MediumInterface {
    // TODO: implement.
    const Medium *inside;
    const Medium *outside;
};

#endif // PBRT_CORE_MEDIUM_H