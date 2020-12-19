#include "sampler.h"

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

Float PixelSampler::Get2D() {
    if (current2DDimension < samples2D.size()) {
        return samples2D[current2DDimension++][currentPixelSampleIndex];
    } else {
        // If the caller keeps requesting samples beyond the number of known dimensions
        // (the nSampledDimensions number that was passed to the constructor and that
        // was used to allocated dimension arrays in samples2D), the sampler returns
        // uniform random numbers.
        return rng.UniformFloat();
    }
}