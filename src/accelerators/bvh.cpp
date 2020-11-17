#include "bvh.h"

struct BVHPrimitiveInfo {
    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;

    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber),
          bounds(bounds),
          centroid(.5f*bounds.pMin + .5f*bounds.pMax) 
    {}
};

struct BVHBuildNode {
    // Bounding box that contains all the primitives in nodes beneath it.
    Bounds3f bounds;

    // Primitives get assigned to one child or the other according to which side of the
    // split axis they are on.
    BVHBuildNode *children[2];

    int splitAxis;

    // The primitives contained by this node are stored contiguously in the
    // BVHAccel::primitives array as a result of the tree's construction.
    int firstPrimOffset;
    int nPrimitives;

    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
    }

    // Axis is the split axis, which is used to decide which child to assign a given
    // primitive.
    void InitInterior(int axis, BVHBuildNode *child0, BVHBuildNode *child1) {
        children[0] = child0;
        children[1] = child1;
        bounds = Union(child0->bounds, child1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
};

BVHAccel(
    const std::vector<std::shared_ptr<Primitive>> &p,
    int maxPrimsInNode = 1,
    SplitMethod splitMethod = SplitMethod::SAH
) : primitives(p), maxPrimsInNode(std::min(255, maxPrimsInNode), splitMethod(splitMethod) {
    if (primitives.size() == 0) {
        return;
    }
    
    // Compute the bounding box (and centroid) of each primitive and collect them.
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
        primitiveInfo[i] = { i, primitives[i]->WorldBound() };
    }

    // Build BVH tree.

    int totalNodes = 0;
    MemoryArena arena(1024 * 1024);

    // Primitives in leaf nodes will occupy contiguous ranges here.
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
    
    BVHBuildNode *root;
    if (splitMethod == SplitMethod::HLBVH) {
        root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
    } else {
        root = recursiveBuild(
            arena, primitiveInfo, 0, primitives.size(), &totalNodes, orderedPrims
        );
    }

    primitives.swap(orderedPrims);

    // TODO: compute representation of depth-first traversal of BVH tree.
}

BVHBuildNode *BVHAccel::recursiveBuild(
    MemoryArena &arena,
    std::vector<BVHPrimitiveInfo> &primitiveInfo,
    int start,
    int end,
    int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims
) {
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();
    (*totalNodes)++;

    // Compute bounds of all primitives in node.
    Bounds3f bounds;
    for (int i = start; i < end; ++i) {
        bounds = Union(bounds, primitiveInfo[i].bounds);
    }

    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf node and store pointers to associated primitives in contiguous range.
        int firstPrimOffset = orderedPrims.size();
        for (int i = start; i < end; ++i) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        
        return node;
    } else {
        // Compute bounds of primitive centroids.
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i) {
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
        }

        // Split primitives into 2 partitions. One partition goes to one child and the
        // other to the other. Ideally, the bounds of one partition doesn't overlap with
        // the other's, but that's not always possible.
        //
        // Partitioning along the (centroid) bound's axis of largest extent often leaves
        // some room between the bounds of the 2 partitions. This is not always true, though:
        // one of the other axes of smaller extent might result in less overlap.
        int dim = centroidBounds.MaximumExtent();

        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Zero-volume centroid bounds. All the primitives' centroids are aligned; note
            // that nPrimitives > 0 here.

            // Create leaf node.
            int firstPrimOffset = orderedPrims.size();
            for (int i = start; i < end; ++i) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            
            return node;
        } else {
            // Partition primitives based on splitMethod. Note that if splitMethod is HLBVH
            // recursiveBuild wouldn't have been invoked by the BVHAccel constructor. 
            switch (splitMethod){
                case SplitMethod::Middle: {
                    // TODO.
                }
                case SplitMethod::EqualCounts: {
                    // TODO.
                }
                case SplitMethod::SAH:
                default: {

                } 
            }

            node->InitInterior(dim,
                recursiveBuild(arena, primitiveInfo, start, mid, totalNodes, orderedPrims),
                recursiveBuild(arena, primitiveInfo, mid, end, totalNodes, orderedPrims)
            );
        }
    }
    
    return node;
}