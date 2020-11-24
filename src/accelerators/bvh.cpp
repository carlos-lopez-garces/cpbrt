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

                // Partition primitives through node's midpoint.
                case SplitMethod::Middle: {
                    // The midpoint of the interval between the centroids that are farthest apart
                    // along the dominant axis.
                    Float pmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                    
                    // Partition the primitiveInfo array into 2 subsets: primitives whose centroids
                    // lie to the left of the midpoint will be put in the 1st partition; in the 2nd one
                    // otherwise. midPtr will mark the end of the 1st partition and the beginning of the
                    // 2nd one.
                    BVHPrimitiveInfo *midPtr = std::partition(
                        &primitiveInfo[start],
                        &primitiveInfo[end-1]+1,
                        [dim, pmid](const BVHPrimitiveInfo &pi) {
                            return pi.centroid[dim] < pmid;
                        }
                    );

                    mid = midPtr - &primitiveInfo[0];
                    if (mid != start && mid != end) {
                        break;
                    }
                    
                    // No partitioning ocurred; this is supposed to occur when all the primitives have
                    // overlapping bounding boxes. Fall through to next split method.
                }

                // Partition primitives into equally sized subsets.
                case SplitMethod::EqualCounts: {
                    // Re-establish the original mid primitive position in case execution fell through from
                    // the previous split method.
                    mid = (start + end) / 2;

                    // Partition the primitiveInfo array into 2 subsets of equal size: primitives whose
                    // centroids lie to the left of the centroid of the primitive in the middle will be
                    // put in the 1st partition; in the 2nd one otherwise. Note that this isn't a full
                    // O(nlogn) sort: there's no order relation among the primitives of a partition, other
                    // than being to the left or to the right of the centroid of the mid primitive. O(n).
                    std::nth_element(
                        &primitiveInfo[start],
                        &primitiveInfo[mid],
                        &primitiveInfo[end-1]+1,
                        [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
                            return a.centroid[dim] < b.centroid[dim];
                        }
                    );
                }

                // Partition primitives using approximate surface area heuristic.
                case SplitMethod::SAH:
                default: {
                    if (nPrimitives <= 4) {
                        // TODO: Partition primitives into equally sized subsets.
                    } else {
                        // Split the dominant axis into a constant number of buckets of equal length.
                        constexpr int nBuckets = 12;
                        struct BucketInfo{
                            int count = 0;
                            Bounds3f bounds;
                        };
                        BucketInfo buckets[nBuckets];

                        // Place primitives into buckets.
                        for (int i = start; i < end; ++i) {
                            // Bounds.Offset returns the position of the primitive's centroid relative
                            // to the extent of the bound's dominant axis. For example, a primitive centroid
                            // lying at 1 quarter of the extent of the bound's x-axis has an offset of 0.25.
                            //
                            // Note that b is an integer such that 0 <= b <= 12 that maps a continuous
                            // offset to one of the 12 discrete buckets.
                            int b = nBuckets * centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
                            
                            if (b == nBuckets) {
                                b = nBuckets - 1;
                            }
                            buckets[b].count++;
                            buckets[b].bounds = Union(buckets[b].bounds, primitiveInfo[i].bounds);
                        }

                        // Stores the cost of each candidate partitioning. 
                        Float cost[nBuckets - 1];

                        // Instead of computing the actual cost of computation, a guess of their relative costs
                        // is used. Ray-primitive intersection is guessed to be 8 times higher than node traversal,
                        // which is basically a ray-AABB intersection test. A single ray-primitive intersection's
                        // cost is given by the primitive's surface area.
                        Float nodeTraversalRelativeCost = 0.125f;
                        Float primitiveIntersectionRelativeCost = 1.0f;

                        // Each candidate partitioning splits the primitives at one of the bucket boundaries.
                        for (int i = 0; i < nBuckets - 1; ++i) {
                            // b0 is the AABB of the 1st partition, b1 is the AABB of the 2nd partition.
                            Bounds3f b0, b1;
                            int count0 = 0, count1 = 0;

                            // 1st partition includes buckets 0 to i.
                            for (int j = 0; j <= i; ++j) {
                                b0 = Union(b0, buckets[j].bounds);
                                count0 += buckets[j].count;
                            }

                            // 2nd partition includes buckets i to nBuckets-1.
                            for (int j = i+1; j < nBuckets; ++j) {
                                b1 = Union(b1, buckets[j].bounds);
                                count1 += buckets[j].count;
                            }

                            // Let 'bounds' be this interior node and let its 2 children be b0 and b1. The conditional
                            // probability that a ray will intersect b0 given that 'bounds' has been intersected is the
                            // ratio of their surface areas. Likewise for b1.
                            //
                            // p(b0|bounds) = b0.SurfaceArea() / bounds.SurfaceArea()
                            // 
                            // p(b1|bounds) = b1.SurfaceArea() / bounds.SurfaceArea()
                            //
                            // So, the ratio approaches 1 as b0.SurfaceArea() approaches bounds.SurfaceArea() (likewise
                            // for b1). The larger the volume of 'bounds' occupied by b1 is, the greater the conditional
                            // probability is.
                            cost[i] 
                                = nodeTraversalRelativeCost
                                  + (
                                      count0 * primitiveIntersectionRelativeCost * b0.SurfaceArea() 
                                      + count1 * primitiveIntersectionRelativeCost * b1.SurfaceArea()
                                  ) / bounds.SurfaceArea();
                        }

                        // Find bucket (partition boundary) to split at that minimizes SAH.
                        Float minCost = cost[0];
                        int minCostSplitBucket = 0;
                        for (int i = 1; i < nBuckets - 1; ++i) {
                            if (cost[i] < minCost) {
                                minCost = cost[i];
                                minCostSplitBucket = i;
                            }
                        }

                        // Either split primitives at selected SAH bucket and continue the recursion, or create leaf node.
                        Float leafCost = nPrimitives * primitiveIntersectionRelativeCost;
                        if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
                            // Even if creating a leaf node has an estimated lesser cost, maxPrimsInNode restricts the
                            // number of primitives, so splitting (at the minimal-cost bucket boundary) is required.
                            BVHPrimitiveInfo *midPtr = std::partition(
                                &primitiveInfo[start],
                                &primitiveInfo[end-1]+1,
                                [=](const BVHPrimitiveInfo &pi) {
                                    int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                    if (b == nBuckets) {
                                        b = nBuckets - 1;
                                    }
                                    return b <= minCostSplitBucket;
                                }
                            );

                            mid = midPtr - &primitiveInfo[0];
                        } else {
                            // Create leaf node.
                            int firstPrimOffset = orderedPrims.size();
                            for (int i = start; i < end; ++i) {
                                int primNum = primitiveInfo[i].primitiveNumber;
                                orderedPrims.push_back(primitives[primNum]);
                            }
                            
                            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                            
                            return node;
                        }
                    }
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