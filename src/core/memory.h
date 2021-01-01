// Allocates a block of the given size aligned to a L1 cache line boundary.
void *AllocAligned(size_t size);

// Allocates the given count of blocks of type T contiguously.
template <typename T> T *AllocAligned(size_t count);

void FreeAligned(void * ptr);