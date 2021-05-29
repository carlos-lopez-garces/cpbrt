#include "sampler.h"
#include "sampling.h"
#include "camera.h"

void Sampler::StartPixel(const Point2i &p) {
    currentPixel = p;
    currentPixelSampleIndex = 0;
    
    // Reset array offsets for next pixel sample.
    array1DOffset = 0;
    array2DOffset = 0;
}

bool Sampler::StartNextSample() {
    // Reset array offsets for next pixel sample.
    array1DOffset = 0;
    array2DOffset = 0;

    return ++currentPixelSampleIndex < samplesPerPixel;
}

bool Sampler::SetSampleNumber(int64_t sampleIndex) {
    // Reset array offsets for next pixel sample.
    array1DOffset = 0;
    array2DOffset = 0;

    currentPixelSampleIndex = sampleIndex;
    return currentPixelSampleIndex < samplesPerPixel;
}

CameraSample Sampler::GetCameraSample(const Point2i &pRaster) {
    CameraSample cs;
    cs.pFilm = ((Point2f) pRaster) + Get2D();
    cs.time = Get1D();
    cs.pLens = Get2D();
    return cs;
}

void Sampler::Request1DArray(int n) {
    samples1DArraySizes.push_back(n);
    sampleArray1D.push_back(std::vector<Float>(n * samplesPerPixel));
}

void Sampler::Request2DArray(int n) {
    samples2DArraySizes.push_back(n);
    sampleArray2D.push_back(std::vector<Point2f>(n * samplesPerPixel));
}

const Float *Sampler::Get1DArray(int n) {
    if (array1DOffset == sampleArray1D.size()) {
        return nullptr;
    }

    return &sampleArray1D[array1DOffset++][currentPixelSampleIndex * n];
}

const Point2f *Sampler::Get2DArray(int n) {
    if (array2DOffset == sampleArray2D.size()) {
        return nullptr;
    }

    return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
}

bool PixelSampler::StartNextSample() {
    current1DDimension = 0;
    current2DDimension = 0;

    return Sampler::StartNextSample();
}

bool PixelSampler::SetSampleNumber(int64_t sampleIndex) {
    current1DDimension = 0;
    current2DDimension = 0;

    return Sampler::SetSampleNumber(sampleIndex);
}

Float PixelSampler::Get1D() {
    if (current1DDimension < samples1D.size()) {
        return samples1D[current1DDimension++][currentPixelSampleIndex];
    } else {
        // If the caller keeps requesting samples beyond the number of known dimensions
        // (the nSampledDimensions number that was passed to the constructor and that
        // was used to allocated dimension arrays in samples1D), the sampler returns
        // uniform random numbers.
        return rng.UniformFloat();
    }
}

Point2f PixelSampler::Get2D() {
    if (current2DDimension < samples2D.size()) {
        return samples2D[current2DDimension++][currentPixelSampleIndex];
    } else {
        // If the caller keeps requesting samples beyond the number of known dimensions
        // (the nSampledDimensions number that was passed to the constructor and that
        // was used to allocated dimension arrays in samples2D), the sampler returns
        // uniform random numbers.
        return Point2f(rng.UniformFloat(), rng.UniformFloat());
    }
}

void GlobalSampler::StartPixel(const Point2i &p) {
    Sampler::StartPixel(p);
    
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(0);

    // Compute arrayEndDim for dimensions used for array samples.
    arrayEndDim = arrayStartDim + sampleArray1D.size() + 2*sampleArray2D.size();
    
    // Compute 1D array samples for GlobalSampler.
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        // samples1DArraySizes records the number of dimensions per sample of a given
        // array request.
        int nSamples = samples1DArraySizes[i] * samplesPerPixel;

        for (int j = 0; j < nSamples; ++j) {
            int64_t index = GetIndexForSample(j);
            sampleArray1D[i][j] = SampleDimension(index, arrayStartDim + i);
        }
    }

    // Compute 2D array samples for GlobalSampler.
    int dim = arrayStartDim + samples1DArraySizes.size();
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        // samples2DArraySizes records the number of dimensions per sample of a given
        // array request.
        int nSamples = samples2DArraySizes[i] * samplesPerPixel;

        for (int j = 0; j < nSamples; ++j) {
            int64_t index = GetIndexForSample(j);
            sampleArray2D[i][j].x = SampleDimension(index, dim);
            sampleArray2D[i][j].y = SampleDimension(index, dim + 1);
        }
    }
}

bool GlobalSampler::StartNextSample() {
    // There's a distinction to make between samples and pixel samples: a sample is the value
    // generated for a given dimension (or 2 at once, if the sample is 2D); a pixel sample is
    // the collection of values generated for all of its dimensions.
    //
    // StartNextSample() starts a new pixel sample. That's why the recorded next dimension is
    // reset.
    dimension = 0;

    intervalSampleIndex = GetIndexForSample(currentPixelSampleIndex + 1);

    return Sampler::StartNextSample();
}

bool GlobalSampler::SetSampleNumber(int64_t sampleNum) {
    // There's a distinction to make between samples and pixel samples: a sample is the value
    // generated for a given dimension (or 2 at once, if the sample is 2D); a pixel sample is
    // the collection of values generated for all of its dimensions.
    //
    // SetSampleNumber() sets the current pixel sample. That's why the recorded next dimension
    // is reset.
    dimension = 0;

    intervalSampleIndex = GetIndexForSample(sampleNum);

    return Sampler::SetSampleNumber(sampleNum);
}

Float GlobalSampler::Get1D() {
    if (dimension >= arrayStartDim && dimension < arrayEndDim) {
        dimension = arrayEndDim;
    }

    return SampleDimension(intervalSampleIndex, dimension++);
}

Point2f GlobalSampler::Get2D() {
    if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim) {
        dimension = arrayEndDim;
    }

    Point2f p(
        SampleDimension(intervalSampleIndex, dimension),
        SampleDimension(intervalSampleIndex, dimension+1)
    );

    dimension += 2;

    return p;
}