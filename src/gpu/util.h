#ifndef CPBRT_GPU_UTIL_H
#define CPBRT_GPU_UTIL_H

#include "core/cpbrt.h"
#include "core/error.h"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>

#define CUDA_CHECK(EXPR)                                    \
    if (EXPR != cudaSuccess) {                              \
        cudaError_t error = cudaGetLastError();             \
        Error("CUDA error: %s", cudaGetErrorString(error)); \
    } else

#endif // CPBRT_GPU_UTIL_H