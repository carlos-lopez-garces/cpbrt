#include <algorithm>
#include <utility>

#include "spectrum.h"

// Determines whether the samples are sorted in order of increasing wavelength.
bool SpectrumSamplesSorted(const Float *lambda, int n) {
    for (int i = 0; i < n-1; ++i) {
        if (lambda[i] > lambda[i+1]) {
            return false;
        }
    }
    return true;
}


// Sorts the samples in order of increasing wavelength.
void SortSpectrumSamples(Float *lambda, Float *values, int n) {
    std::vector<std::pair<Float, Float>> sortArray;
    sortArray.reserve(n);
    for (int i = 0; i < n; ++i) {
        sortArray.push_back(std::make_pair(lambda[i], values[i]));
    }
    std::sort(sortArray.begin(), sortArray.end());
    for (int i = 0; i < n; ++i) {
        lambda[i] = sortArray[i].first;
        values[i] = sortArray[i].second;
    }
}

// Computes the average value of the samples contained in the wavelength subinterval
// [lambdaStart, lambdaEnd]. This subinterval is one of the nSpectralSamples subintervals of
// the regular partition of the entire domain [sampledLambdaStart, sampledLambdaEnd] of 
// wavelengths.
//
// If the SPD function was a known continuous one, its average value over the [lambdaStart, lambdaEnd]
// subinterval would be the area under the curve divided by the length of the subinterval. The area
// under the curve would be given by the definite integral of the function on the [lambdaStart, lambdaEnd]
// interval of integration.
//
// But we only have a sampled discretization of the function from which we have to reconstruct it.
// The reconstruction is a piecewise linear function of line segments between samples. The average over
// the subinterval is still the area under the curve divided by the length of the subinterval, but the
// area under the curve can only be approximated by the sum of the areas of rectangles between samples
// (like a general Riemman sum over a general partition, but with possibly multiple rectangles per
// subinterval). 
Float AverageSpectrumSamples(
    const Float *lambda, 
    const Float *values,
    int n,
    Float lambdaStart,
    Float lambdaEnd
) {
    // Handle cases with out-of-bounds range or single sample only.

    // The subinterval precedes the first sample, [start, end, 0th sample].
    // The SPD function is assumed to be constant on the subinterval and equal to the first sample.
    if (lambdaEnd <= lambda[0]) return values[0];

    // The subinterval follows the last sample, [n-1th sample, start, end].
    // The SPD function is assumed to be constant on the subinterval and equal to the last sample.
    if (lambdaStart >= lambda[n-1]) return values[n-1];

    // Regardless of the sample's wavelength, if the wavelength interval is not partitioned, the
    // SPD function is regarded to be constant on the entire domain [sampledLambdaStart, sampledLambdaEnd]
    // and equal to the only sample's value.
    if (n == 1) return values[0];

    Float sum = 0;
    
    // Add contributions of constant segments before/after samples.

    // [start, 0th sample, end], between the start of the subinterval and the first sample, the SPD
    // function is regarded to be constant and equal to the first sample's value. Its contribution
    // to the sum is equal to the area under the curve in the subinterval [start, 0th sample]; since
    // the function is constant here, the region under the curve is rectangular.
    if (lambdaStart < lambda[0]) {
        sum += values[0] * (lambda[0] - lambdaStart);
    }

    // [start, end, n-1th sample], between the last sample and the end of the subinterval, the SPD
    // function is regarded to be constant and equal to the last sample's value. Its contribution
    // to the sum is equal to the area under the curve in the subinterval [n-1th sample, end]; since
    // the function is constant here, the region under the curve is rectangular.
    if (lambdaEnd > lambda[n-1]) {
        sum += values[n-1] * (lambdaEnd - lambda[n-1]);
    }
    
    // A segment is the line segment that connects 2 samples. These segments are the linear
    // reconstruction of the sampled SPD function.
    //
    // Find the 1st segment that straddles the left endpoint of the [lambdaStart, lambdaEnd] subinterval,
    // if any: [lambda[i], lambdaStart, lambda[i+1], lambdaEnd].
    int i = 0;
    while (lambdaStart > lambda[i+1]) {
        ++i;
    }

    // The lerp lambda function interpolates linearly between the values of 2 adjacent samples
    // that form a segment. A segment may lie completely inside the subinterval [lambdaStart, lambdaEnd]
    // or it may straddle one of the endpoints:
    // [lambda[i], start, lambda[i+1], end] or [start, lambda[i], end, lambda[i+1]].
    //
    // When the segment lies inside the subinterval, w will be equal to either lambda[i] or lambda[i+1]
    // and, in those cases, no interpolation is needed because the values of the SPD function at those
    // wavelengths are known (the sample).
    //
    // When the segment straddles 1 or the 2 endpoints of the subinterval, w is lambda[i] clamped to
    // lambdaStart and/or lambda[i+1] clampled to lambdaEnd. In those cases, the value of the SPD 
    // function is not known and needs to be obtained by linear interpolation.
    auto lerp = [lambda, values](Float w, int i) {
        // When w is equal to lambda[i], the interpolator is 0 and lambda[i] is itself the output of
        // the interpolation. Likewise for lambda[i+1], but there the interpolator is 1.
        return Lerp(
            (w - lambda[i]) / (lambda[i+1] - lambda[i]),
            values[i],
            values[i+1]
        );
    };

    // Loop over wavelength sample segments and add contributions. 
    for (; i+1 < n && lambdaEnd >= lambda[i]; ++i) {
        // One of the segments may straddle the left endpoint of the [lambdaStart, lambdaEnd]
        // subinterval; in that case, lambda[i] is clamped to lambdaStart.
        //
        // Similarly, one of the segments may straddle the right endpoint; in that case, lambda[i+1]
        // is clamped to lambdaEnd.
        //
        // In yet another case, a single segment may straddle the 2 endpoints of the subinterval,
        // [lambda[i], lambdaStart, lambdaEnd, lambda[i+1]]; in that case, lambda[i] is also clamped
        // to lambdaStart and lambda[i+1] is clamped to lambdaEnd.
        Float segmentLambdaStart = std::max(lambdaStart, lambda[i]);
        Float segmentLambdaEnd = std::min(lambdaEnd, lambda[i+1]);

        // The contribution of this segment is the midpoint value of the values of the 2 (possibly
        // clampled) samples divided by the length of the (possibly clampled) distance between the 2.
        // This amount corresponds to the area of the rectangle that exists between the 2 samples
        // (its height being the midpoint value of the value of the 2 samples). 
        sum += 0.5 
            * (lerp(segmentLambdaStart, i) + lerp(segmentLambdaEnd, i))
            * (segmentLambdaEnd - segmentLambdaStart);
    }

    return sum / (lambdaEnd - lambdaStart);
}