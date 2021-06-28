#ifndef CPBRT_ACCELERATORS_BVH_H
#define CPBRT_ACCELERATORS_BVH_H

#include <atomic>
#include <memory>
#include <vector>

#include "core/cpbrt.h"
#include "core/primitive.h"

struct BVHBuildNode;
struct BVHPrimitiveInfo;
struct LinearBVHNode;
struct MortonPrimitive;

class BVHAccel : public Aggregate {
public:
    // Primitive subdivision algorithms.
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    BVHAccel(
        const std::vector<std::shared_ptr<Primitive>> &p,
        int maxPrimsInNode = 1,
        SplitMethod splitMethod = SplitMethod::SAH
    );

    Bounds3f WorldBound() const;

    BVHBuildNode *recursiveBuild(
        MemoryArena &arena,
        // Array of unprocessed primitive bounding boxes and centroids.
        std::vector<BVHPrimitiveInfo> &primitiveInfo,
        // The range [start, end) of the primitiveInfo array that this recursive call handles.
        // The result is a partitioning of the primitives into [start, mid) and [mid, end),
        // which branch out the recursion further.
        int start,
        int end,
        int *totalNodes,
        // Pointers to primitives processed by this subtree of the recursion will be 
        // stored contiguously in the corresponding range of this array. This is actually how
        // an interior node of the BVH refers to its primitives: with an index range over this
        // array.
        std::vector<std::shared_ptr<Primitive>> &orderedPrims
    );

    BVHBuildNode *HLBVHBuild(
        MemoryArena &arena,
        const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims
    ) const;

    BVHBuildNode *emitLBVH(
        BVHBuildNode *&buildNodes,
        const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        MortonPrimitive *mortonPrims,
        int nPrimitives,
        int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims,
        std::atomic<int> *orderedPrimsOffset,
        int bitIndex
    ) const;

    BVHBuildNode *buildUpperSAH(
        MemoryArena &arena,
        std::vector<BVHBuildNode *> &treeletRoots,
        int start,
        int end,
        int *totalNodes
    ) const;

    int flattenBVHTree(BVHBuildNode *node, int *offset);

    bool Intersect(const Ray &ray, SurfaceInteraction *si) const;

    bool IntersectP(const Ray &ray) const;

private:
    // Maximum number of primitives that can be stored in a node, up to 255.
    const int maxPrimsInNode;

    // Algorithm to use for subdividing primitives.
    const SplitMethod splitMethod;

    // Before construction, the primitives lie in an undefined order. After construction,
    // the order is determined by the split method. The primitives of a leaf node will
    // lie contiguously.
    std::vector<std::shared_ptr<Primitive>> primitives;

    // Binary tree in linear array representation and depth-first, built during construction. 
    LinearBVHNode *nodes = nullptr;
};

std::shared_ptr<BVHAccel> CreateBVHAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

#endif // CPBRT_ACCELERATORS_BVH_H