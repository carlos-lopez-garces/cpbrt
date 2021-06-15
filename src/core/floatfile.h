#ifndef PBRT_CORE_FLOATFILE_H
#define PBRT_CORE_FLOATFILE_H

#include "cpbrt.h"

bool ReadFloatFile(const char *filename, std::vector<Float> *values);

#endif // PBRT_CORE_FLOATFILE_H