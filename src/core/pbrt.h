#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Choose between float and double for most floating-point declarations
// at compile time.
#ifdef CPBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // CPBRT_FLOAT_AS_DOUBLE