#include <vector>

#include "cpbrt.h"
#include "geometry.h"
#include "rng.h"

void StratifiedSample1D(Float *sample, int nSamples, RNG &rng, bool jitter);

void StratifiedSample2D(Point2f *sample, int xSamples, int ySamples, RNG &rng, bool jitter);

// Randomly permutes an array of samples.
template <typename T> void Shuffle(T *sample, int nSamples, int nDimensions, RNG &rng) {
    for (int i = 0; i < nSamples; ++i) {
        int other = i + rng.UniformUInt32(nSamples - i);

        for (int j = 0; j < nDimensions; ++j) {
            std::swap(
                sample[nDimensions * i + j],
                sample[nDimensions * other + j]
            );
        }
    }
}

void LatinHypercube(Float *sample, int nSamples, int nDimensions, RNG &rng);

struct Distribution1D {
    // The n images of a piecewise-constant function.
    std::vector<Float> f;
    
    // The definite integral of f on [0,1].
    Float definiteIntegral;

    // f's cumulative distribution function CDF corresponding to the probability distribution
    // function PDF obtained from f (f itself is not a PDF, because its image values may be
    // unbounded and its definite integral over the domain may not be 1).
    std::vector<Float> cdf;

    Distribution1D(
        // n constant values of a piecewise-constant function f with domain [0,1].
        const Float *f,
        int n
    ) : f(f, f + n), cdf(n+1) {
        cdf[0] = 0;
        // Compute definite integral of step function at xi.
        //
        // i is the 0-based index of the ith subinterval of the domain [0,1], i=0,1,...,n.
        // Note that this doesn't compute the CDF as such, but the definite integral of f
        // at the overlapping subintervals [0,i=1], [0,i=2], ..., [0,i=n]: for a piecewise-
        // constant function, or step function, like f, the definite integral is the sum of
        // the products of the step's value and the (constant) step length: SUM(f(i)*1/n).
        //
        // The definite integral of f on its entire domain [0,1] will end up in cdf[n].
        //
        // This is just a step towards computing the CDF, a step that happens to compute the
        // definite integral of f already.  
        for (int i = 1; i < n+1; ++i) {
            cdf[i] = cdf[i-1] + f[i-1] / n;
        }

        definiteIntegral = cdf[n];

        if (definiteIntegral == 0) {
            for (int i = 1; i < n+1; ++i) {
                cdf[i] = Float(i) / Float(n);
            }
        } else {
            // Transform integral of step function into CDF.
            //
            // It can be shown that the PDF p(x) for a piecewise-constant function like f 
            // with definite integral F is p(x) = f(x)/F. Intuitively, think of the probability
            // of xi as the fraction of the total area under the curve of f that the ith 
            // rectangle f(xi)*deltaX represents (as deltaX->0, as in a Riemann sum).
            //
            // The CDF P(x) is obtained by integrating the PDF. By dividing by the definite integral
            // of f, we are implicitly turning the intermediary cdf[i] value into the sum of the
            // PDF evaluated at xi, p(xi), and the CDF at xi-1; and explicitly turning it into
            // the CDF evaluated at xi: P(xi) = P(xi-1) + p(xi).
            //
            // The resulting CDF is a continuous piecewise-linear function, where each subinterval
            // is a linear function that starts where the last segment left off, and of slope f(xi)/nF:
            // P(xi) = P(xi-1) + p(xi) = P(xi-1) + f(xi)/nF. 
            // TODO: why the n?
            for (int i = 1; i < n+1; ++i) {
                cdf[i] /= definiteIntegral;
            }
        }
    }

    // Samples the distribution determined by f, returning a sample x and its probability p(x) as assigned
    // by the corresponding PDF. In essence, the sample x is found by finding the inverse of the CDF P and
    // evaluating it at a [0,1] uniform random sample u. The draw of x will obey the probability distribution
    // determined by f. This is known as the Inverse Transform Method of generating random variables of a
    // desired distribution.
    Float SampleContinuous(
        // [0,1] uniform random sample to be used to sample this PDF. 
        Float u,
        Float *pdf,
        int *off = nullptr
    ) const {
        // Each index i of cdf[], i=0,1...n, is a grid point of a regular partition of the [0,1] domain
        // that results in n subintervals [xi-1, xi]. cdf[i] maps the ith subinterval [xi-1, xi] to the
        // cumulative probability of xi, the right endpoint of the subinterval: cdf[i] = P(xi).
        //
        // We want to map the [0,1] uniform sample u to the domain value x that has cumulative probability
        // u: x such that P(X <= x) = u. Since cdf[] only records the cumulative probability of the endpoints,
        // the point x with P(X <= x) = u must be found somewhere in the subinterval [xi, xi+1], xi <= x < xi+1,
        // where cdf[i] <= u. The offset found next is that index i.  
        int i = FindInterval(
            cdf.size(),
            [&](int index) {
                return cdf[index] <= u;
            }
        );

        if (off) *off = i;

        // Compute sampled offset.
        //
        // The point x with cumulative probability P(X <= x) = u lies in [xi, xi+1): cdf[i] <= u < cdf[i+1].
        // The exact cumulative probability of x is found by linearly interpolating cdf[i] and cdf[i+1], 
        // because this CDF is piecewise-linear. Because of this linearity, the same resulting interpolated
        // offset corresponds to x.
        //
        // Here, the cumulative probability of x is cdf[i] + du. And x = xi + du.
        Float du = u - cdf[i];
        if ((cdf[i+1] - cdf[i]) > 0) {
            du /= cdf[i+1] - cdf[i];
        }
        Float dx = du;

        // Obtain and evaluate PDF for sampled x: p(x).
        if (pdf) {
            *pdf = f[i] / definiteIntegral;   
        }

        // Return x. This was all equivalent to finding the inverse of the CDF and evaluating it at u.
        return (i + dx) / Count();
    }

    // Samples the discrete distribution determined by f, which assigns a probability to the ith step
    // (f is piecewise-constant, or a step function) that is proportional to f(i) and equal to the 
    // probability mass function PMF determined by f.
    int SampleDiscrete(Float u, Float *pdf = nullptr, Float *uRemapped = nullptr) const {
        int i = FindInterval(
            cdf.size(),
            [&](int index) {
                return cdf[index] <= u;
            }
        );

        if (pdf) {
            *pdf = f[i] / (definiteIntegral * Count());
        }

        if (uRemapped) {
            *uRemapped = (u - cdf[i]) / (cdf[i+1] - cdf[i]);
        }

        return i;
    }

    int Count() const {
        return f.size();
    }
};