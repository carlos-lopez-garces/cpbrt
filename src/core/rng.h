#include "cpbrt.h"

// Largest representable floating-point number that is less than 1.
#ifdef CPBRT_FLOAT_AS_DOUBLE
static const Float OneMinusEpsilon = 0x1.fffffffffffffp-1;
#else
static const Float OneMinusEpsilon = 0x1.fffffep-1;
#endif