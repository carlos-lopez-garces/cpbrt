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

void *MemoryArena::Alloc(size_t nBytes) {
    // Round up requested size to multiple of 16 bytes. A 16-byte-aligned pointer
    // ensures that any primitive data type will be word-aligned, including the
    // largest one, long double (16 bytes in 64-bit architectures). 
    nBytes = ((nBytes + 15) & (~15));

    if (currentBlockPos + nBytes > currentAllocSize) {
        // No available memory in currentBlock to satisfy request.
        
        // Add current block to usedBlocks list.
        if (currentBlock) {
            usedBlocks.push_back(std::make_pair(currentAllocSize, currentBlock));
            currentBlock = nullptr;
        }

        // Try to get block from availableBlocks with enough memory to satisfy request.
        for (auto iter = availableBlocks.begin(); iter != availableBlocks.end(); ++iter) {
            if (iter->first >= nBytes) {
                currentAllocSize = iter->first;
                currentBlock = iter->second;
                availableBlocks.erase(iter);
                break;
            }
        }

        if (!currentBlock) {
            // No available block to reuse. Allocate a new one.
            currentAllocSize = std::max(nBytes, blockSize);
            currentBlock = AllocAligned<uint8_t>(currentAllocSize);
        }
        currentBlockPos = 0;
    }

    void *requested = currentBlock + currentBlockPos;
    currentBlockPos += nBytes;
    return requested;
}

template <typename T> T *MemoryArena::Alloc(size_t n, bool runConstructor) {
    T *requested = (T *) Alloc(sizeof(T) * n);

    if (runConstructor) {
        for (size_t i = 0; i < n; ++i) {
            new (&requested[i]) T();
        }
    }

    return requested;
}