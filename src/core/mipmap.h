#ifndef CPBRT_CORE_MIPMAP_H
#define CPBRT_CORE_MIPMAP_H

#include "cpbrt.h"

// What to do when the supplied texture coordinates are outside [0,1].
enum class ImageWrap {
    Repeat,
    Black,
    Clamp
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

            // Resample in the s direction.

            // Resample in the t direction.
        }
    }
};

#endif // CPBRT_CORE_MIPMAP_H
