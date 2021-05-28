#ifndef CPBRT_CORE_MEMORY_H
#define CPBRT_CORE_MEMORY_H

#include <cstdint>
#include <list>
#include <utility>

#include "cpbrt.h"

// Allocates an object of the given type using the supplied MemoryArena.
// The use of the new operator is for calling the constructor on the memory
// allocated by the MemoryArena, not to allocate it.
#define ARENA_ALLOC(arena, Type) new (arena.Alloc(sizeof(Type))) Type

// Allocates a block of the given size aligned to a L1 cache line boundary.
void *AllocAligned(size_t size);

// Allocates the given count of blocks of type T contiguously.
template <typename T> T *AllocAligned(size_t count);

void FreeAligned(void * ptr);

class MemoryArena {
private:
    // Minimum unit of allocation of real memory.
    const size_t blockSize;
 
    // Block where allocations are made from currently.
    uint8_t *currentBlock = nullptr;

    // Offset of first free location in currentBlock.
    size_t currentBlockPos = 0;

    // Size of currentBlock (at least blockSize, but may be larger if an allocation
    // request is made for >blockSize bytes).
    size_t currentAllocSize = 0;

    // Past blocks that are now full.
    std::list<std::pair<size_t, uint8_t *>> usedBlocks;

    // Reset blocks that are now available.
    std::list<std::pair<size_t, uint8_t *>> availableBlocks;

public:
    MemoryArena(size_t blockSize = 262144 /* 256kB*/) : blockSize(blockSize) {}

    // Allocates a 16-byte-aligned object of the given size in the currentBlock. 
    void *Alloc(size_t nBytes);

    // Allocates the given count of objects of type T contiguously.
    template <typename T> T *Alloc(size_t n = 1, bool runConstructor = true);

    // Reset the current and all the past blocks.
    void Reset() {
        currentBlock = 0;
        availableBlocks.splice(availableBlocks.begin(), usedBlocks);
    }
};

#endif // CPBRT_CORE_MEMORY_H