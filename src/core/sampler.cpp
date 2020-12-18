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