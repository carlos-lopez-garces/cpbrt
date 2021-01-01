#include <malloc.h>

#include "cpbrt.h"
#include "memory.h"

// Allocates a block of the given size aligned to a L1 cache line boundary.
void *AllocAligned(size_t size) {
#if defined(CPBRT_HAVE_ALIGNED_MALLOC)
    return _aligned_malloc(size, CPBRT_L1_CACHE_LINE_SIZE);
#else
    return malloc(size);
#endif
}

// Allocates the given count of blocks of type T contiguously.
template <typename T> T *AllocAligned(size_t count) {
    return (T *) AllocAligned(sizeof(T) * count);
}

void FreeAligned(void *ptr) {
    if (!ptr) {
        return;
    }

#if defined(CPBRT_HAVE_ALIGNED_MALLOC)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}