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
template <typename T> T *AllocAligned(size_t count) {
    return (T *) AllocAligned(sizeof(T) * count);
}

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
    void *Alloc(size_t nBytes) {
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
                currentAllocSize = 0;
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

    // Allocates the given count of objects of type T contiguously.
    template <typename T> T *Alloc(size_t n = 1, bool runConstructor = true) {
        T *requested = (T *) Alloc(sizeof(T) * n);

        if (runConstructor) {
            for (size_t i = 0; i < n; ++i) {
                new (&requested[i]) T();
            }
        }

        return requested;
    }

    // Reset the current and all the past blocks.
    void Reset() {
        currentBlockPos = 0;
        availableBlocks.splice(availableBlocks.begin(), usedBlocks);
    }
};

// Blocked arrays subdivide the data into square blocks of a small fixed power-of-2 size.
// Entries are identified with a (u,v) coordinate. Each u coordinate maps to a block and
// then an offset into that block; the v coordinate too: (b_u, b_v) is the block number 
// of a given (u,v) coordinate and (o_u, o_v) is the offset into it. These block numbers
// and offsets are used to compute the index of the (u,v) into the linear data array.
template <typename T, int logBlockSize> class BlockedArray {
private:
    T *data;
    const int uRes;
    const int vRes;
    const int uBlocks;

public:
    // Entries in the 
    // logBlockSize is the exponent of the power of 2 that is the size of each block.
    BlockedArray(int uRes, int vRes, const T *d = nullptr)
      : uRes(uRes), vRes(vRes), uBlocks(RoundUp(uRes) >> logBlockSize) {

        int nAlloc = RoundUp(uRes) * RoundUp(vRes);
        data = AllocAligned<T>(nAlloc);

        for (int i = 0; i < nAlloc; ++i) {
            // "Placement" syntax: specify the location where the object of type T
            // will be initialized.
            new (&data[i]) T(); 
        }

        if (d) {
            for (int v = 0; v < vRes; ++v) {
                for (int u = 0; u < uRes; ++u) {
                    // Copy from the source data into the corresponding block.
                    (*this)(u, v) = d[v * uRes + u];
                }
            }
        }
    }

    constexpr int BlockSize() const {
        // 2^logBlockSize.
        return 1 << logBlockSize;
    }

    int RoundUp(int x) const {
        return (x + BlockSize() - 1) & ~(BlockSize() - 1);
    }

    int uSize() const {
        return uRes;
    }

    int vSize() const {
        return vRes;
    }

    // (u,v) entries are indexed by block first and then offset.
    int Block(int a) const {
        return a >> logBlockSize;
    }

    int Offset(int a) const {
        return (a & (BlockSize() - 1));
    }

    // Retrieves the entry at the (u,v) coordinate from the linear array.
    T &operator()(int u, int v) {
        // Compute block number.
        int bu = Block(u), bv = Block(v);
        // Compute offset into block.
        int ou = Offset(u), ov = Offset(v);
        // Compute offset into linear array.
        int offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
        offset += BlockSize() * ov + ou;
        return data[offset];
    }
};

#endif // CPBRT_CORE_MEMORY_H