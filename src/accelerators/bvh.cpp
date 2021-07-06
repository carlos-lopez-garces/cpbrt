#include "bvh.h"
#include "core/interaction.h"
#include "core/paramset.h"
#include "core/parallel.h"

struct BVHPrimitiveInfo {
    // Index of primitive in BVHAccel::primitives.
    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;

    BVHPrimitiveInfo() {}

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

struct LinearBVHNode {
    // AABB that encompasses this subtree's or leaf's primitives.
    Bounds3f bounds;

    union {
        // If node is leaf.
        int primitivesOffset;
        // If node is interior, its first child lies right next to it, so there's no
        // need to store its offset, only the offset to its second child.
        int secondChildOffset;
    };

    // Interior nodes don't refer to primitives directly, only to their 2 children
    // nodes. Only leaf nodes refer to primitives directly. 
    uint16_t nPrimitives;
    
    // The axis along which the primitives were partitioned at this node, x, y, or z.
    // Only for interior nodes.
    uint8_t axis;

    // Make efficient use of cache by ensuring the struct is 32 bytes in size: if this
    // node is cache-line aligned, the next one in the array will also be, and none
    // will straddle cache lines.
    uint8_t pad[1]; 
};

struct LBVHTreelet {
    // Linear array index where the primitives of this Morton cluster / treelet start.
    int startIndex;
    int nPrimitives;
    BVHBuildNode *buildNodes;
};

struct MortonPrimitive {
    // Index of the primitive in the primitiveInfo array.
    int primitiveIndex;
    uint32_t mortonCode;
};

// Takes a 32-bit value and returns the result of shifting the ith bit to be at the 
// 3ith bit, leaving zeros in other bits.
inline uint32_t LeftShift3(uint32_t x) {
    // TODO: define.
    // CHECK_LE(x, (1 << 10));
    if (x == (1 << 10)) --x;
#ifdef CPBRT_HAVE_BINARY_CONSTANTS
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0b00000011000000001111000000001111;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0b00000011000011000011000011000011;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0b00001001001001001001001001001001;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
    x = (x | (x << 16)) & 0x30000ff;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0x300f00f;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0x30c30c3;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0x9249249;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // CPBRT_HAVE_BINARY_CONSTANTS
    return x;
}

// The Morton encoding of a coordinate interleaves the bits of the z, y, and x components,
// ...z3y3x3z2y2x2z1y1x1.
inline uint32_t EncodeMorton3(const Vector3f &v) {
    return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
}

static void RadixSort(std::vector<MortonPrimitive> *mortonPrims) {
    std::vector<MortonPrimitive> tmpMortonPrims(mortonPrims->size());
    constexpr int bitsPerPass = 6;
    constexpr int nBits = 30;
    constexpr int nPasses = nBits / bitsPerPass;

    // Sort bitsPerPass bits per pass. A 30-bit long Morton code is sorted digit by digit, starting
    // at the low bits, 6 bits at a time. It takes 5 passes to sort the 30 bits 6 bits at a time.
    //
    // To sort in passes, 2 arrays are needed: mortonPrims and tmpMortonPrims; in the 1st pass,
    // mortonPrims is the unsorted input array and tmpMortonPrims will be a 6-bit-sorted
    // output array; in the 2nd pass, tmpMortonPrims is the partially-sorted input array and
    // mortonPrims will be a 12-bit-sorted output array; etc. The alternation continues pass after
    // pass.
    for (int pass = 0; pass < nPasses; ++pass) {
        int lowBit = pass * bitsPerPass;

        // 'in' is a reference to the input array to be sorted. It alternates on every pass.
        std::vector<MortonPrimitive> &in = (pass & 1) ? tmpMortonPrims : *mortonPrims;
        // 'out' is a reference to the array that will store the sorted values. It alternates too.
        std::vector<MortonPrimitive> &out = (pass & 1) ? *mortonPrims : tmpMortonPrims;

        constexpr int bitMask = (1 << bitsPerPass) - 1;

        // 2^bitsPerPass buckets.
        constexpr int nBuckets = 1 << bitsPerPass;
        // Count the Morton primitives that will land on each bucket.
        int bucketCount[nBuckets] = { 0 };
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            ++bucketCount[bucket];
        }

        // The output array will store the Morton primitives of a given bucket contiguously.
        // The bucket counters are used to calculate the starting index of 'out' at which
        // the Morton primitives of a given bucket will be stored.
        int outIndex[nBuckets];
        outIndex[0] = 0;
        for (int i = 1; i < nBuckets; ++i) {
            outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];
        }

        // Store sorted Morton primitives at the calculated indices of the output array.
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            // As Morton primitives of a given bucket get stored at the bucket's output index,
            // the bucket's output index needs to get advanced for the next primitive of the
            // same bucket.
            out[outIndex[bucket]++] = mp;
        }
    }

    // Due to the alternation, the fully-sorted array may be mortonPrims or tmpMortonPrims.
    if (nPasses & 1) {
        std::swap(*mortonPrims, tmpMortonPrims);
    }
}

BVHAccel::BVHAccel(
    const std::vector<std::shared_ptr<Primitive>> &p,
    int maxPrimsInNode,
    SplitMethod splitMethod
) : primitives(p), maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod) {
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

    // Primitives in leaf nodes will occupy contiguous intervals here.
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

    // Build linear array representation of binary tree in depth-first order.
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    int offset = 0;
    flattenBVHTree(root, &offset);
}

Bounds3f BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : Bounds3f();
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
                        [dim](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) {
                            return a.centroid[dim] < b.centroid[dim];
                        }
                    );
                }

                // Partition primitives using approximate surface area heuristic.
                case SplitMethod::SAH:
                default: {
                    if (nPrimitives <= 2) {
                        // Partition the primitiveInfo array into 2 subsets of equal size: primitives whose
                        // centroids lie to the left of the centroid of the primitive in the middle will be
                        // put in the 1st partition; in the 2nd one otherwise. Note that this isn't a full
                        // O(nlogn) sort: there's no order relation among the primitives of a partition, other
                        // than being to the left or to the right of the centroid of the mid primitive. O(n).
                        mid = (start + end) / 2;
                        std::nth_element(
                            &primitiveInfo[start],
                            &primitiveInfo[mid],
                            &primitiveInfo[end-1]+1,
                            [dim](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) {
                                return a.centroid[dim] < b.centroid[dim];
                            }
                        );
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
                            // to the extent of the bound's along each of the axes. For example, a primitive
                            // centroid lying at 1 quarter of the extent of the bound's x-axis has an offset of 0.25.
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

BVHBuildNode *BVHAccel::HLBVHBuild(
    MemoryArena &arena,
    const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims
) const {
    // Compute bounding box of all primitive centroids.
    Bounds3f bounds;
    for (const BVHPrimitiveInfo &pi : primitiveInfo) {
        bounds = Union(bounds, pi.centroid);
    }

    // Compute Morton codes of primitives. Each worker thread is given 512 primitives to
    // process. 
    std::vector<MortonPrimitive> mortonPrims(primitiveInfo.size());
    ParallelFor(
        [&](int i) {
            // Bottom 10 bits of x, bottom 10 bits of y, and bottom 10 bits of z so that the total
            // 30 bits fit in a 32-bit variable.
            constexpr int mortonBits = 10;

            // 2^10.
            constexpr int mortonScale = 1 << mortonBits;
            
            mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;

            Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid);

            // Bounds offsets are floating-point numbers in [0, 1], but the Morton transformation takes 
            // a coordinate with integer components, so the offsets are discretized by mapping them
            // to regularly-spaced integers in the range [0, 1024] by scaling them by mortonScale=2^10.
            mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
        },
        primitiveInfo.size(),
        512
    );

    // Radix sort Morton primitives. Sorted Morton codes have special properties that allow for spatial
    // reasoning about the relative positions of the primitives.
    RadixSort(&mortonPrims);

    // Create treelets for Morton clusters (primitives with nearby centroids along the Morton curve)
    // at bottom of BVH.
    //
    // Take a binary string of 3 bits as example. Let the 2 upper bits be 2^2  groups: 00X, 01X, 10X, and 11X.
    // Each group has members XX0 and XX1.
    //
    // Now take the 30-bit Morton code binary strings. Let the 12 upper bits be 2^12 groups. The remaining
    // low 18 bits represent the 2^18 members of each of the 2^12 groups.
    //
    // Spatially, a group is a cell of a grid in 3-space, a subvolume.
    //
    // All the primitives of a Morton cluster lie in the same cell.

    // Find intervals of primitives for each treelet.
    std::vector<LBVHTreelet> treeletsToBuild;
    for (int start = 0, end = 1; end <= (int)mortonPrims.size(); ++end) {
#ifdef CPBRT_HAVE_BINARY_CONSTANTS
        uint32_t mask = 0b00111111111111000000000000000000;
#else
        uint32_t mask = 0x3ffc0000;
#endif

        // Morton codes of the same cell have the same 12 upper bits. A cell is wrapped up when the 12 upper bits
        // change from the last Morton primitive to the current one, or when the current one is the last one. A
        // treelet is created for the cluster in the interval [start, end] of mortonPrims.
        if (end == (int)mortonPrims.size()
            || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask))
        ) {
            // Allocate treelet and add it to the array.
            int nPrimitives = end - start;
            int maxBVHNodes = 2 * nPrimitives;

            // Passing false to the allocator prevents the constructor of BVHBuildNode from being invoked. 
            BVHBuildNode *nodes = arena.Alloc<BVHBuildNode>(maxBVHNodes, false);
            // For the time being, the treelet's buildNodes pointer will point at the start of the interval of
            // the cluster on mortonPrims's memory. This pointer will later point to the root of the LBVH that
            // will be built for this treelet.
            treeletsToBuild.push_back({start, nPrimitives, nodes});

            start = end;
        }
    }

    // Create LBVHs for treelets in parallel. 
    //
    // atomicTotal and orderedPrimsOffset will be updated by parallel or concurrent threads, so they need to be
    // accessed atomically.
    //
    // The lower [17, 0] bits split the primitives at an ever finer level of granularity [17 -> 0]. At these levels,
    // the primitives are always split at the midpoint of the axis that corresponds to the Morton code bit
    // (...ZYXZYX...ZYX).
    std::atomic<int> atomicTotal(0);
    std::atomic<int> orderedPrimsOffset(0);
    orderedPrims.resize(primitives.size());
    ParallelFor(
        [&](int i) {
            // Generate ith LBVH treelet.
            int nodesCreated = 0;
            const int firstBitIndex = 29 - 12;
            LBVHTreelet &treelet = treeletsToBuild[i];
            treelet.buildNodes = emitLBVH(
                treelet.buildNodes,
                primitiveInfo,
                &mortonPrims[treelet.startIndex],
                treelet.nPrimitives,
                &nodesCreated,
                orderedPrims,
                &orderedPrimsOffset,
                firstBitIndex
            );
            atomicTotal += nodesCreated;
        },
        treeletsToBuild.size()
    );
    *totalNodes = atomicTotal;

    std::vector<BVHBuildNode *> finishedTreelets;
    for (LBVHTreelet &treelet : treeletsToBuild) {
        finishedTreelets.push_back(treelet.buildNodes);
    }

    // The upper [29, 16] bits split at an ever lower level of granularity [29 <- 16]. At these levels,
    // treelets are split using the surface area heuristic (SAH).
    return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(), totalNodes);
}

// emitLBVH is invoked on a Morton cluster. All the primitives lie in the same grid cell or
// subvolume in 3-space.
BVHBuildNode *BVHAccel::emitLBVH(
    BVHBuildNode *&buildNodes,
    const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    MortonPrimitive *mortonPrims,
    int nPrimitives,
    int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims,
    // Initially 0.
    std::atomic<int> *orderedPrimsOffset,
    // Initially 17, the first low bit that is not a group.
    int bitIndex
) const {
    // bitIndex counts down from 17 to -1. Every time it decreases, the primitives get further
    // subdivided by 2 and split into 2 groups:
    //
    // bitIndex=17: XXXXXXXXXXXX0................. and XXXXXXXXXXXX1.................
    // bitIndex=16: XXXXXXXXXXXXX0................ and XXXXXXXXXXXXX1................
    // ...
    // bitIndex=0:  XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0 and XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1
    //
    // Since a Morton code interleaves the bits of the coordinate's components ZYXZYX...ZYX,
    // a different split axis is used every time bitIndex decreases. Primitives where the bit's
    // value is 0 go in one group and the ones with 1 go in the other group. If all of the
    // primitives have the same value at the bit, the primitives won't be split along this axis
    // at this level, and the algorithm proceeds to the next lower bit.
    if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {
        // Create and return leaf node for LBVH treelet of this split Morton subcluster.
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++;
        Bounds3f bounds;
        // fetch_add atomically reads the current offset and adds the input count to it.
        // This effectively reserves an interval of orderedPrims addresses for the current
        // Morton subcluster while avoiding data races with other concurrent or parallel
        // threads.
        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
        for (int i = 0; i < nPrimitives; ++i) {
            int primitiveIndex = mortonPrims[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex];
            bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        int mask = 1 << bitIndex;

        // Advance to next lower bit if all the primitives lie on the same side of this bit's axis.
        if ((mortonPrims[0].mortonCode & mask) == (mortonPrims[nPrimitives-1].mortonCode & mask)) {
            return emitLBVH(
                buildNodes,
                primitiveInfo,
                mortonPrims,
                nPrimitives,
                totalNodes,
                orderedPrims,
                orderedPrimsOffset,
                bitIndex - 1
            );
        }
        
        // Find split point for this bit's axis using a binary search. The split point is the unique, but
        // possibly non-existent, point where this bit's value changes from 0 to 1 across Morton codes.
        // All the primitives to the left of the split point have a bit value of 0, whereas the ones to
        // the right have a 1.
        int searchStart = 0;
        int searchEnd = nPrimitives-1;
        while (searchStart+1 != searchEnd) {
            int mid = (searchStart + searchEnd) / 2;
            if ((mortonPrims[searchStart].mortonCode & mask) == (mortonPrims[mid].mortonCode & mask)) {
                searchStart = mid;
            } else {
                searchEnd = mid;
            }
        }
        int splitOffset = searchEnd;

        // Create interior LBVH node that splits along this bit's axis the primitives into its 2 children.
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++;
        BVHBuildNode *lbvh[2] = {
            emitLBVH(
                buildNodes,
                primitiveInfo,
                mortonPrims,
                splitOffset,
                totalNodes,
                orderedPrims,
                orderedPrimsOffset,
                bitIndex-1
            ),
            emitLBVH(
                buildNodes,
                primitiveInfo,
                &mortonPrims[splitOffset],
                nPrimitives-splitOffset,
                totalNodes,
                orderedPrims,
                orderedPrimsOffset,
                bitIndex-1
            )
        };
        int axis = bitIndex % 3;
        node->InitInterior(axis, lbvh[0], lbvh[1]);
        return node;
    }
}

BVHBuildNode *BVHAccel::buildUpperSAH(
    MemoryArena &arena,
    std::vector<BVHBuildNode *> &treeletRoots,
    int start, int end,
    int *totalNodes
) const {
    // TODO: define.
    // CHECK_LT(start, end);
    int nNodes = end - start;
    if (nNodes == 1) {
        return treeletRoots[start];
    }
    (*totalNodes)++;
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();

    // Compute bounds of all nodes under this HLBVH node.
    Bounds3f bounds;
    for (int i = start; i < end; ++i) {
        bounds = Union(bounds, treeletRoots[i]->bounds);
    }

    // Compute bound of HLBVH node centroids, choose split dimension _dim_.
    Bounds3f centroidBounds;
    for (int i = start; i < end; ++i) {
        Point3f centroid = (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) * 0.5f;
        centroidBounds = Union(centroidBounds, centroid);
    }
    int dim = centroidBounds.MaximumExtent();
    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something?
    // TODO: define.
    // CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

    // Allocate _BucketInfo_ for SAH partition buckets.
    CPBRT_CONSTEXPR int nBuckets = 12;
    struct BucketInfo {
        int count = 0;
        Bounds3f bounds;
    };
    BucketInfo buckets[nBuckets];

    // Initialize _BucketInfo_ for HLBVH SAH partition buckets.
    for (int i = start; i < end; ++i) {
        Float centroid = (treeletRoots[i]->bounds.pMin[dim] + treeletRoots[i]->bounds.pMax[dim]) * 0.5f;
        int b = nBuckets * ((centroid - centroidBounds.pMin[dim]) / (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
        if (b == nBuckets) {
            b = nBuckets - 1;
        }
        // TODO: define.
        // CHECK_GE(b, 0);
        // CHECK_LT(b, nBuckets);
        buckets[b].count++;
        buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
    }

    // Compute costs for splitting after each bucket.
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
        Bounds3f b0, b1;
        int count0 = 0, count1 = 0;
        for (int j = 0; j <= i; ++j) {
            b0 = Union(b0, buckets[j].bounds);
            count0 += buckets[j].count;
        }
        for (int j = i + 1; j < nBuckets; ++j) {
            b1 = Union(b1, buckets[j].bounds);
            count1 += buckets[j].count;
        }
        cost[i] = .125f + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
    }

    // Find bucket to split at that minimizes SAH metric.
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }

    // Split nodes and create interior HLBVH SAH node.
    BVHBuildNode **pmid = std::partition(
        &treeletRoots[start], &treeletRoots[end - 1] + 1,
        [=](const BVHBuildNode *node) {
            Float centroid = (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
            int b = 
                nBuckets * ((centroid - centroidBounds.pMin[dim]) / (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
            if (b == nBuckets) b = nBuckets - 1;
            // TODO: define.
            // CHECK_GE(b, 0);
            // CHECK_LT(b, nBuckets);
            return b <= minCostSplitBucket;
        }
    );
    int mid = pmid - &treeletRoots[0];
    // TODO: define.
    // CHECK_GT(mid, start);
    // CHECK_LT(mid, end);
    node->InitInterior(
        dim, this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
        this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes)
    );
    return node;
}

// Takes a linked representation of the BVH binary tree and flattens it into a linear
// array binary tree in depth-first order.
int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];

    // This effectively places this node at this offset of the linear array binary tree.
    linearNode->bounds = node->bounds;

    // Remember this node's offset (to return it) and advance it for the next recursive
    // calls.
    int nodeOffset = (*offset)++;

    if (node->nPrimitives > 0) {
        // Leaf node.
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Interior node.
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        
        // The first child subtree will lie on a contiguous segment just next to this
        // node's offset.
        //
        // Upon return, offset will be the next one after this child subtree's last
        // linear node.
        flattenBVHTree(node->children[0], offset);

        // The second child subtree will lie after all of the linear nodes of the first
        // child subtree. The returned offset (that of the second child node itself) is
        // recorded in this linear node.
        linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
    }

    return nodeOffset;
}

// Traverses the BVH to find the closest intersection.
bool BVHAccel::Intersect(const Ray &ray, SurfaceInteraction *si) const {
    bool hit = false;
    Vector3f reciprocalDir(1/ray.d.x, 1/ray.d.y, 1/ray.d.z);
    int dirIsNeg[3] = { reciprocalDir.x < 0, reciprocalDir.y < 0, reciprocalDir.z < 0 };

    // toVisitOffset points at the top of nodesToVisit, which acts as a stack (of offsets
    // to linear nodes). When traversing an interior node, the closer of the 2 children will
    // be visited immediately and the farther one will be pushed onto the stack to be
    // processed next as soon as the closer child has been traversed completely.
    int toVisitOffset = 0;
    int nodesToVisit[64];
    int currentNodeIndex = 0;

    // Follow ray through BVH (linear) nodes to find primitive intersections.
    while (true) {
        const LinearBVHNode *node = &nodes[currentNodeIndex];

        // TODO: implement Bounds3::IntersectP overload.
        if (node->bounds.IntersectP(ray) /*node->bounds.IntersectP(ray, reciprocalDir, dirIsNeg)*/) {
            // The ray intersects the node's AABB.

            if (node->nPrimitives > 0) {
                // Leaf node. Test primitives for intersection (recall that only leaf nodes
                // store primitives).
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // The ray's extent, tMax, is shortened every time a closer intersection is
                    // found. Subsequent farther intersections get efficiently discarded.
                    if (primitives[node->primitivesOffset + i]->Intersect(ray, si)) {
                        hit = true;
                    }
                    
                    // Keep going, this intersection might not be the closest one.
                }

                if (toVisitOffset == 0) {
                    break;
                }

                // Next node to visit.
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Interior node. First, determine which of its 2 children the ray intersects
                // first (if at all). By traversing the closer child first, it is more likely that
                // the closest intersection will be found, which would shorten the ray's extent,
                // tMax, causing subsequent intersections with the farther child to be discarded
                // efficiently.
                // 
                // Put farther BVH node on nodesToVisit stack.
                if (dirIsNeg[node->axis]) {
                    // Visit the second child node first. Put the other one on the stack.
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    // Visit the first child node first. Put the other one on the stack.
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            // The ray didn't intersect the node's AABB.

            if (toVisitOffset == 0) {
                break;
            }

            // Next node to visit.
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    return hit;
}

// Tells whether the ray intersects a primitive. The traversal of the BVH stops at the
// first primitive intersection.
bool BVHAccel::IntersectP(const Ray &ray) const {
    Vector3f reciprocalDir(1/ray.d.x, 1/ray.d.y, 1/ray.d.z);
    int dirIsNeg[3] = { reciprocalDir.x < 0, reciprocalDir.y < 0, reciprocalDir.z < 0 };

    // toVisitOffset points at the top of nodesToVisit, which acts as a stack (of offsets
    // to linear nodes). When traversing an interior node, the closer of the 2 children will
    // be visited immediately and the farther one will be pushed onto the stack to be
    // processed next as soon as the closer child has been traversed completely without having
    // found an intersection.
    int toVisitOffset = 0;
    int nodesToVisit[64];
    int currentNodeIndex = 0;

    // Follow ray through BVH (linear) nodes to find the 1st primitive intersection.
    while (true) {
        const LinearBVHNode *node = &nodes[currentNodeIndex];

        if (node->bounds.IntersectP(ray, reciprocalDir, dirIsNeg)) {
            // The ray intersects the node's AABB.

            if (node->nPrimitives > 0) {
                // Leaf node. Test primitives for intersection (recall that only leaf nodes
                // store primitives).
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // One intersection is enough.
                    if (primitives[node->primitivesOffset + i]->IntersectP(ray)) {
                        return true;
                    }
                }

                if (toVisitOffset == 0) {
                    break;
                }

                // Next node to visit.
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Interior node. First, determine which of its 2 children the ray intersects
                // first (if at all). By traversing the closer child first, it is more likely that
                // the first intersection will be found quicker.
                // 
                // Put farther BVH node on nodesToVisit stack.
                if (dirIsNeg[node->axis]) {
                    // Visit the second child node first. Put the other one on the stack.
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    // Visit the first child node first. Put the other one on the stack.
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            // The ray didn't intersect the node's AABB.

            if (toVisitOffset == 0) {
                break;
            }

            // Next node to visit.
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    return false;
}

// TODO: explain.
std::shared_ptr<BVHAccel> CreateBVHAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
    BVHAccel::SplitMethod splitMethod;
    if (splitMethodName == "sah")
        splitMethod = BVHAccel::SplitMethod::SAH;
    else if (splitMethodName == "hlbvh")
        splitMethod = BVHAccel::SplitMethod::HLBVH;
    else if (splitMethodName == "middle")
        splitMethod = BVHAccel::SplitMethod::Middle;
    else if (splitMethodName == "equal")
        splitMethod = BVHAccel::SplitMethod::EqualCounts;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".", splitMethodName.c_str());
        splitMethod = BVHAccel::SplitMethod::SAH;
    }

    int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return std::make_shared<BVHAccel>(std::move(prims), maxPrimsInNode, splitMethod);
}