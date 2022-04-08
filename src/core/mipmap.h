#ifndef CPBRT_CORE_MIPMAP_H
#define CPBRT_CORE_MIPMAP_H

#include "cpbrt.h"
#include "texture.h"
#include "parallel.h"

// What to do when the supplied texture coordinates are outside [0,1].
enum class ImageWrap {
    Repeat,
    Black,
    Clamp
};

// Return type of MIPMap::resampleWeights.
struct ResampleWeight {
    // Offset to the first of the 4 original texels that contribute to the value
    // of the given texel in the resampled, higher-resolution image. These 4 texels
    // are contiguous.
    int firstTexel;
    // 4 texels of the original image contribute to the value of 1 pixel in a
    // resampled, higher-resolution image.
    Float weight[4];
};

// The structure of the MIPMap is a pyramid where each level, from the base up, has
// increasingly lower resolution. Each texel of each level is the result of filtering
// 4 texels of the lower, higher-resolution level. The base level is the original image
// (or the resampled one if its resolution wasn't a power of 2).
template <typename T> class MIPMap {
private:
    const bool doTrilinearFiltering;
    const Float maxAnisotropy;
    const ImageWrap wrapMode;
    std::vector<std::unique_ptr<BlockedArray<T>>> pyramid;

public:
    template <typename T> MIPMap(
        const Point2i &resolution,
        const T *image,
        bool doTrilinearFiltering,
        Float maxAnisotropy,
        ImageWrap wrapMode
    ) : doTrilinearFiltering(doTrilinearFiltering),
        maxAnisotropy(maxAnisotropy),
        wrapMode(wrapMode),
        resolution(resolution) {

        std::unique_ptr<T[]> resampledImage = nullptr;
        if (!IsPowerOf2(resolution.x) || !IsPowerOf2(resolution.y)) {
            // Ensure that the resolution of the image is a power of 2 in each
            // dimension (s and t), by resampling it.
            Point2i resolutionPow2(RoundUpPow2(resolution.x), RoundUpPow2(resolution.y));

            // We are going to resample an image of resolution (s, t) at a higher
            // sampling rate to obtain an image of resolution (s^2, t^2). We'll use
            // a separable filter f(s,t)=f(s)f(t), so  we can use a 1D filter to filter
            // in the s direction and obtain a resized image of resolution (s^2, t) and
            // then resample the result in the t direction to obtain the final image of
            // resolution (s^2, t^2). 

            // s-pass: Resample (magnify) in the s direction. The resulting image will have a
            // resolution of (resolutionPow2[0], resolution[1]).
            std::unique_ptr<ResampleWeight[]> sWeights = resampleWeights(resolution.x, resolutionPow2.x);
            // Allocate space for the final image of resolution (resolutionPow2[0], resolutionPow2[1]),
            // not just for this pass in the s direction.
            resampledImage.reset(new T[resolutionPow2.x * resolutionPow2.y]);
            ParallelFor(
                // Each t is one row of texels of the original image.
                [&](int t) {
                    // Each s is a texel in this row.
                    for (int s = 0; i < resolutionPow2.x; ++s) {
                        // 0 corresponds to ImageWrap::Black. It'll be overwritten if wrapMode is
                        // actually Repeat or Clamp.
                        resampledImage[t * resolutionPow2.x + s] = 0.f;
                        // Each texel of the resampled, higher-resolution image takes
                        // contributions from 4 texels of the original image, the ones that
                        // surround the resampled sample point in the original image (the center
                        // of the reconstruction filter).
                        for (int j = 0; j < 4; ++j) {
                            int originalS = sWeights[s].firstTexel + j;

                            if (wrapMode == ImageWrap::Repeat) {
                                originalS = Mod(originalS, resolution.x);
                            } else if (wrapMode == ImageWrap::Clamp) {
                                originalS = Clamp(originalS, 0, resolution.x-1);
                            }

                            if (originalS >= 0 && originalS < (int) resolution.x) {
                                // New texel in the resampled, higher-resolution image is the weighted
                                // average of 4 original texels that surround the image look up point. 
                                resampledImage[t * resolutionPow2.x + s] += sWeights[s].weight[j] * image[t * resolution.x + originalS];
                            }
                        }
                    }
                },
                // Each thread processes 16 of the texels in the given row.
                resolution.y,
                16
            );

            // t-pass: Resample (magnify) in the t direction.
            std::unique_ptr<ResampleWeight[]> tWeights = resampleWeights(resolution.y, resolutionPow2.y);
            std::vector<T *> resampleBuffers;
            int nThreads = MaxThreadIndex();
            // Give each thread a column of the resized image to work on.
            for (int i = 0; i < nThreads; ++i) {
                resampleBuffers.push_back(new T[resolutionPow2.y]);
            }
            ParallelFor(
                // Each s is one column of texels of the original image.
                [&](int s) {
                    T *workData = resampleBuffers[ThreadIndex];

                    // Each t is a texel in this column.
                    for (int t = 0; t < resolutionPow2.y; ++t) {
                        // 0 corresponds to ImageWrap::Black. It'll be overwritten if wrapMode is
                        // actually Repeat or Clamp.
                        workData[t] = 0.f;
                        // Each texel of the resampled, higher-resolution image takes
                        // contributions from 4 texels of the original image, the ones that
                        // surround the resampled sample point in the original image (the center
                        // of the reconstruction filter).
                        for (int j = 0; j < 4; ++j) {
                            int offset = tWeights[t].firstTexel + j;

                            if (wrapMode == ImageWrap::Repeat) {
                                offset = Mod(offset, resolution.y);
                            } else if (wrapMode == ImageWrap::Clamp) {
                                offset = Clamp(offset, 0, (int) resolution.y-1);
                            }

                            if (offset >= 0 && offset < (int) resolution.y) {
                                // New texel in the resampled, higher-resolution image is the weighted
                                // average of 4 original texels that surround the image look up point.
                                // The weight influences the value of the new texel that was set in the
                                // s-pass.
                                workData[t] += tWeights[t].weight[j] * resampledImage[offset * resolutionPow2.x + s];
                            }
                        }
                    }

                    // Set the final texel values.
                    for (int t = 0; t < resolutionPow2.y; ++t) {
                        resampledImage[t * resolutionPow2.x + s] = Clamp(workData[t], 0.f, Infinity);
                    }
                },
                // Each thread processes 32 of the texels in the given column.
                resolutionPow2.x,
                32
            );

            for (auto buffer : resampleBuffers) delete[] buffer;

            resolution = resolutionPow2;
        }
        
        // Initialize mip map levels from image.
        int nLevels = 1 + Log2Int(std::max(resolution.x, resolution.y));
        pyramid.resize(nLevels);

        // Initialize the most detailed level of the MIPMap (the base), that is, the original
        // image or the resampled one if its resolution wasn't a power of 2.
        pyramid[0].reset(new BlockedArray<T>(resolution.x, resolution.y, resampledImage ? resampledImage.get() : image));

        // Initialize the rest of the levels.
        for (int i = 1; i < nLevels; ++i) {
            int sRes = std::max(1, pyramid[i-1]->uSize() / 2);
            int tRes = std::max(1, pyramid[i-1]->vSize() / 2);
            pyramid[i].reset(new BlockedArray<T>(sRes, tRes));
            // Filter 4 texels from the previous, higher-resolution level to give each texel
            // of this level its value.
            ParallelFor(
                [&](int t) {
                    for (int s = 0; s < sRes; ++s) {
                        (*pyramid[i])(s, t) = .25f * (Texel(i-1, 2*s, 2*t)   + Texel(i-1, 2*s+1, 2*t) + Texel(i-1, 2*s, 2*t+1) + Texel(i-1, 2*s+1, 2*t+1));
                    }
                }, 
                tRes, 
                16
            );
        }

        // TODO: initialize EWA filter weights if needed.
    }

    int Width() const {
        return resolution.x;
    }

    int Height() const {
        return resolution.y;
    }

    int Levels() const {
        return pyramid.size();
    }

    // Returns the value of the texel at (s,t). The ImageWrap choice applies.
    template <typename T> const T &Texel(int level, int s, int t) const;

private:
    // resampleWeights returns one ResampleWeight object per original texel in column or row
    // of the original image. The information returned, which is for one column or one row,
    // applies to the rest of the columns or rows of the image.
    std::unique_ptr<ResampleWeight[]> resampleWeights(int oldResolution, int newResolution) {
        std::unique_ptr<ResampleWeight> weight(new ResampleWeight[newResolution]);
        // 2 original texels. 4 texels end up contributing to each new texel in the resampled,
        // higher resolution image: 2 in the s and 2 in the t direction.
        Float filterWidth = 2.f;
        for (int i = 0; i < newResolution; ++i) {
            // Center of the filter in the original image. Since the filter's width is 2,
            // the 4 texels that surround this center contribute to the value of the new
            // texel in the resampled, higher resolution image.
            Float center = (i+0.5f) * oldResolution / newResolution;
            // Index of the 1st of the 4 original texels that contribute to this new texel.
            weight[i].firstTexel = std::floor((center - filterWidth) + 0.5f);
            for (int j = 0; j < 4; ++j) {
                Float originalTexelPosition = weight[i].firstTexel + j + 0.5f;
                weight[i].weight[j] = LanczosFilter((originalTexelPosition - center) / filterWidth);
            }

            // Normalize.
            Float reciprocalSumWeights = 1 / (weight[i].weight[0] + weight[i].weight[1] + weight[i].weight[2] + weight[i].weight[3]);
            for (int j = 0; j < 4; ++j) {
                // Normalize the weight of each contributing texel, because the sum might not
                // sum to 1 (depends on the filter).
                weight[i].weight[j] *= reciprocalSumWeights;
            }
        }
        return weight;
    }
};

#endif // CPBRT_CORE_MIPMAP_H
