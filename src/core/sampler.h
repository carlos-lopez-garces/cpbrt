#ifndef CPBRT_CORE_SAMPLER_H
#define CPBRT_CORE_SAMPLER_H

#include <cstdint>
#include <memory>
#include <vector>

#include "cpbrt.h"
#include "geometry.h"
#include "rng.h"

// Samples are collections of multi-dimensional vectors of at least 5 dimensions, (x, y, t, u, v, ...):
// image/film/raster coordinates x and y, time, and lens coordinates u and v.
//
// Samplers aren't aware of the meaning of each dimension. They don't need to: every dimension is a
// [0,1) floating-point number, without exception.
class Sampler {
private:
    // Offset into sampleArray1D of the start of the array of samples of the current dimension. 
    size_t array1DOffset;

    // Offset into sampleArray2D of the start of the array of samples of the current dimension.
    size_t array2DOffset;

protected:
    Point2i currentPixel;

    // From [0, samplesPerPixel). See sampleArray1D and sampleArray2D.
    int64_t currentPixelSampleIndex;

    // Collection of array-samples of 1D samples. An array-sample is a multi-1D-valued sample.
    // The size of each array-sample is recorded in samples1DArraySizes.
    //
    // Sizes: [m, n, p]
    //
    // Array-samples: [
    //      [1D sample 1, ... 1D sample m],
    //      [1D sample 1, ... 1D sample n],
    //      [1D sample 1, ... 1D sample p],
    // ] 
    std::vector<int> samples1DArraySizes;
    std::vector<std::vector<Float>> sampleArray1D;

    // Collection of array-samples of 1D samples. See sampleArray1D.
    std::vector<int> samples2DArraySizes;
    std::vector<std::vector<Point2f>> sampleArray2D;

public:
    const int64_t samplesPerPixel;

    Sampler(int64_t samplesPerPixel) : samplesPerPixel(samplesPerPixel) {}

    // p is the coordinate of the pixel in the image.
    //
    // Base Sampler implementation must be called by subclass implementations.
    virtual void StartPixel(const Point2i &p);

    // Sets the state of the sampler so that it moves on to the first dimension of the next
    // sample for the current pixel at the next request. Returns false after the sampler
    // has been requested to generate n=samplesPerPixel samples; returns true while fewer
    // than n=samplesPerPixel have been requested.
    //
    // Base Sampler implementation must be called by subclass implementations.
    virtual bool StartNextSample();

    // Set the index of the sample to generate next for the current pixel.
    //
    // Base Sampler implementation must be called by subclass implementations.
    virtual bool SetSampleNumber(int64_t sampleIndex);

    // Get the value of the next dimension of the current sample vector.
    virtual Float Get1D() = 0;

    // Get the values of the next 2 dimensions of the current sample vector.
    virtual Point2f Get2D() = 0;

    CameraSample GetCameraSample(const Point2i &pRaster);

    // Allocates an array-sample of size n. Most samplers can generate higher-quality
    // sequences of samples if they generate them at once rather than across a series of
    // Get1D() calls, by accounting for the distribution of sample values across all
    // elements of the array and across the samples in a pixel. 
    void Request1DArray(int n);

    // Allocates an array-sample of size n. Most samplers can generate higher-quality
    // sequences of samples if they generate them at once rather than across a series of
    // Get2D() calls, by accounting for the distribution of sample values across all
    // elements of the array and across the samples in a pixel.
    void Request2DArray(int n);

    // Adjusts the input array size to the sampler's ideal size. Calls to Request1DArray()
    // and Request2DArray() should be made with the size returned by RoundCount(), as that
    // is the number of samples that the particular sampler can generate at once with
    // near-optimal distribution. 
    virtual int RoundCount(int n) const {
        return n;
    }

    // Returns a pointer to the start of the array-sample requested previously with a
    // Request1DArray() call for the same size n.
    const Float *Get1DArray(int n);

    // Returns a pointer to the start of the array-sample requested previously with a
    // Request2DArray() call for the same size n.
    const Point2f *Get2DArray(int n);

    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
};

// PixelSampler generates at once all of the dimensions of all of the vectors of all of the samples
// of a pixel.
//
// PixelSampler delegates the implementation of StartPixel() to its subclasses. These subclasses
// should fill in the samples1D, samples2D, sampleArray1D, and sampleArray2D arrays.
class PixelSampler : public Sampler {
protected:
    int current1DDimension = 0;

    // Collection of arrays of samplesPerPixel 1D samples, organized by dimension.
    //
    // [
    //    [dimension 1: for sample 1, ..., for sample samplesPerPixel], ... [dimension m: for sample 1, ..., for sample samplesPerPixel],
    //    [dimension 1: for sample 1, ..., for sample samplesPerPixel], ... [dimension n: for sample 1, ..., for sample samplesPerPixel],
    //    [dimension 1: for sample 1, ..., for sample samplesPerPixel], ... [dimension p: for sample 1, ..., for sample samplesPerPixel]
    // ]
    //
    // currentPixelSampleIndex is with respect to the array of a given dimension.
    std::vector<std::vector<Float>> samples1D;

    int current2DDimension = 0;

    // Collection of arrays of samplesPerPixel 2D samples, organized by dimension. See samples1D.
    std::vector<std::vector<Point2f>> samples2D;

    RNG rng;

public:
    PixelSampler(int64_t samplesPerPixel, int nSampledDimensions) : Sampler(samplesPerPixel) {
        for (int i = 0; i < nSampledDimensions; ++i) {
            samples1D.push_back(std::vector<Float>(samplesPerPixel));
            samples2D.push_back(std::vector<Point2f>(samplesPerPixel));
        }
    }

    // Overriden.
    bool StartNextSample();
    bool SetSampleNumber(int64_t sampleIndex);
    Float Get1D();
    Point2f Get2D();
};

// GlobalSampler generates consecutive samples that are spread across the entire image, for
// completely different pixels in succession.
class GlobalSampler : public Sampler {
private:
    // Next dimension to generate.
    int dimension;

    // ?
    int64_t intervalSampleIndex;

    // An individual 1D or 2D sample may be requested, or whole arrays of them. Some samplers
    // generate samples of lesser and lesser quality; for example, the sample for the 1st
    // dimension could be of higher quality than that of the 6th.
    //
    // Now, samples for CameraSamples are very important. Thus, the first 5 dimensions are
    // the ones that a Camera uses because they are expected to be high-quality.
    //
    // These dimensions for the Camera are requested by the unit, so the first 5 dimensions
    // are requested in 1s or 2s (via Get1D() and Get2D()). Requests by array are made only
    // over the subsequent dimensions; dimensions [arrayStartDim, arrayEndDim], specifically.
    // Dimensions after arrayEndDim are again for unit requests.
    static const int arrayStartDim = 5;
    int arrayEndDim;

public:
    GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}

    // These 2 functions adapt the base Sampler interface to the needs of a GlobalSampler.
    // A GlobalSampler indexes its samples differently and these functions map between the 2
    // indexing schemes.
    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
    virtual Float SampleDimension(int64_t index, int dimension) const = 0;

    // Overriden.
    void StartPixel(const Point2i &p);
    bool StartNextSample();
    bool SetSampleNumber(int64_t sampleIndex);
    Float Get1D();
    Point2f Get2D();
};

#endif // CPBRT_CORE_SAMPLER_H