#include "stratified.h"
#include "core/sampling.h"

void StratifiedSampler::StartPixel(const Point2i &p) {
    // Generate single stratified 1D samples for the pixel across all dimensions:
    // samples1D[dimension][sample index].
    for (size_t i = 0; i < samples1D.size(); ++i) {
        // Generate all the samples for the current dimension.
        //
        // The number of strata is the same as the total number of pixel samples:
        // xPixelSamples * yPixelSamples, making 1/(xPixelSamples * yPixelSamples)
        // the normalized size of each stratum.
        //
        // The sample index determines the stratum: the kth sample comes from the
        // kth stratum. (Samples will be shuffled next, though.)
        StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng, jitterSamples);
        
        // Randomly permute the samples of the current dimension. Otherwise, there
        // would be a correlation among the values of a given pixel sample across
        // its dimensions: they would all be drawn from the same stratum (from the
        // same [0,1) subinterval). For example, sample kth = (dimension 1: from stratum k,
        // dimension 2: from stratum k, ... dimension n: from stratum k).
        Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }
    
    // Generate single stratified 2D samples for the pixel.
    for (size_t i = 0; i < samples2D.size(); ++i) {
        // Generate all the samples for the current dimension.
        StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng, jitterSamples);
        
        // Randomly permute the samples of the current dimension to avoid correlation.
        Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }

    // Generate arrays of stratified 1D samples for the pixel.
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        int arraySampleSize = samples1DArraySizes[i];
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            // Each array-sample has n=nSamples samples. The 1st array-sample spans indices
            // [0, arraySampleSize); the 2nd array-sample spans [arraySampleSize, 2*arraySampleSize),
            // etc.
            StratifiedSample1D(&sampleArray1D[i][j * arraySampleSize], arraySampleSize, rng, jitterSamples);
            Shuffle(&sampleArray1D[i][j * arraySampleSize], arraySampleSize, 1, rng);
        }
    }

    // Generate arrays of stratified 2D samples for the pixel.
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        int arraySampleSize = samples2DArraySizes[i];
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            LatinHypercube(&sampleArray2D[i][j * arraySampleSize].x, arraySampleSize, 2, rng);
        }
    }

    PixelSampler::StartPixel(p);
}

std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
    StratifiedSampler *ss = new StratifiedSampler(*this);
    ss->rng.SetSequence(seed);
    return std::unique_prt<Sampler>(ss);
}