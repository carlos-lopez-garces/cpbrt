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

template <typename T> class MIPMap {
private:
    const bool doTrilinearFiltering;
    const Float maxAnisotropy;
    const ImageWrap wrapMode;

public:
    template <typename T> MIPMap<T>::MIPMap(
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

            // Resample (magnify) in the s direction. The resulting image will have a
            // resolution of (resolutionPow2[0], resolution[1]).
            std::unique_ptr<ResampleWeight[]> sWeights = resampleWeights(resolution.x, resolutionPow2.x);
            // Allocate space for the final image of resolution (resolutionPow2[0], resolutionPow2[1]),
            // not just for this pass in the s direction.
            resampledImage.reset(new T[resolutionPow2.x * resolutionPow2.y]);

            // Resample (magnify) in the t direction.
        }
    }

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
