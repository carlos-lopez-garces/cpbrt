#include <cstdint>
#include <memory>
#include <vector>

#include "camera.h"
#include "geometry.h"

// Samples are multi-dimensional vectors of at least 5 dimensions, (x, y, t, u, v, ...):
// image/film/raster coordinates x and y, time, and lens coordinates u and v. 
class Sampler {
private:
    // Offset into sampleArray1D of the start of the array of values of the current sample. 
    size_t array1DOffset;

    // Offset into sampleArray2D of the start of the array of values of the current sample.
    size_t array2DOffset;

protected:
    Point2i currentPixel;

    // From [0, samplesPerPixel).
    int64_t currentPixelSampleIndex;

    // Collection of arrays of 1D samples, where a sample is an array of 1D values. The size
    // of each array in sampleArray1D is recorded in samples1DArraySizes and corresponds to the 
    // _number of 1D samples PER sample PER pixel_.
    //
    // Sizes: [m, n, p]
    // Arrays: [
    //    [sample 1: value 1, ..., value m], ... [sample samplesPerPixel: value 1, ..., value m],
    //    [sample 1: value 1, ..., value n], ... [sample samplesPerPixel: value 1, ..., value n],
    //    [sample 1: value 1, ..., value p], ... [sample samplesPerPixel: value 1, ..., value p]
    // ]
    std::vector<int> samples1DArraySizes;
    std::vector<std::vector<Float>> sampleArray1D;

    // Collection of arrays of 2D samples, where a sample is an array of 2D values. The size
    // of each array in sampleArray2D is recorded in samples2DArraySizes and corresponds to the 
    // _number of 2D samples PER sample PER pixel_.
    //
    // Sizes: [m, n, p]
    // Arrays: [[(sample 1), ..., (sample m)], [(sample 1), ..., (sample n)], [(sample 1), ..., (sample p)]]
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

    // Generates an array of n samples per sample per pixel of the next dimension. Most
    // samplers can generate higher-quality sequences of samples if they generate them at
    // once rather than across a series of Get1D() calls, by accounting for the distribution
    // of sample values across all elements of the array and across the samples in a pixel. 
    void Request1DArray(int n);

    // Generates an array of n samples per sample per pixel of the next 2 dimensions. Most
    // samplers can generate higher-quality sequences of samples if they generate them at
    // once rather than across a series of Get2D() calls, by accounting for the distribution
    // of sample values across all elements of the array and across the samples in a pixel.
    void Request2DArray(int n);

    // Adjusts the input array size to the sampler's ideal size. Calls to Request1DArray()
    // and Request2DArray() should be made with the size returned by RoundCount(), as that
    // is the number of samples that the particular sampler can generate at once with
    // near-optimal distribution. 
    virtual int RoundCount(int n) const {
        return n;
    }

    // Returns a pointer to the start of the array requested previously with a Request1DArray()
    // call for the same size n and the current dimension.
    const Float *Get1DArray(int n);

    // Returns a pointer to the start of the array requested previously with a Request2DArray()
    // call for the same size n and the current dimension.
    const Point2f *Get2DArray(int n);

    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
};