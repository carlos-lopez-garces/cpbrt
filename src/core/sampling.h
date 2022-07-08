#ifndef CPBRT_CORE_SAMPLING_H
#define CPBRT_CORE_SAMPLING_H

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

// Samples the unit disk uniformly using a pair of [0,1] uniform random numbers.
// The output sample (X,Y) is the cartesian coordinate of a point on the disk, obtained
// by transforming the polar coordinate of a sample of a 2D random variable (r, theta).
//
// The value of r comes from sampling the marginal PDF of r, which is accomplished by the
// inversion method by evaluating the inverse of the CDF of r with the value of a uniform
// random variable. (The CDF of r is obtained by integrating the marginal PDF.) 
//
// The value of theta comes from sampling the conditional density function p(theta|r) given
// the value of r obtained before. p(theta|r) is sampled by the inversion method by evaluating
// the inverse with the value of a uniform random variable.
Point2f UniformSampleDisk(const Point2f &u);

// Samples the unit disk using a concentric mapping of the unit square, which transforms a
// uniformly distributed random point on the unit square to a point on the unit disk. 
Point2f ConcentricSampleDisk(const Point2f &u);

// TODO: explain.
Vector3f UniformSampleHemisphere(const Point2f &u);

// TODO: explain.
Float UniformHemispherePdf();

// TODO: explain.
Vector3f UniformSampleSphere(const Point2f &u);

// TODO: explain.
Float UniformSpherePdf();

// TODO: explain.
Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax);

// TODO: explain.
Vector3f UniformSampleCone(
    const Point2f &u,
    Float cosThetaMax,
    const Vector3f &x,
    const Vector3f &y,
    const Vector3f &z
);

// TODO: explain.
Float UniformConePdf(Float cosThetaMax);

// Transforms a distribution of points over the unit disk to one of points over the unit
// hemisphere above it, and returns a sample direction.
inline Vector3f CosineSampleHemisphere(const Point2f &u) {
    Point2f d = ConcentricSampleDisk(u);
    Float z = std::sqrt(std::max((Float) 0, 1 - d.x*d.x - d.y*d.y));
    return Vector3f(d.x, d.y, z);
}

// Samples a point of the surface of a triangle uniformly at random, using the input sample
// of a uniform random variable. Sampling is done with a PDF equal to the reciprocal of the
// triangle's area. Note in the implementation, however, that sampling is independent of the
// triangle's dimensions.
//
// The returned point (b0, b1) is a barycentric coordinate that can be used to obtain a
// world-space coordinate on any triangle of any dimension and shape; you'd use it as follows:
//
// p = b0*p0 + b1*p1 + (1 - b0 - b1)*p2
//
// where p0, p1, and p2 are the coordinates of the triangle's vertices. Since barycentric
// coordinates add up to 1, the 3rd component can be derived from the other 2: b2 = 1 - b0 - b1.
Point2f UniformSampleTriangle(const Point2f &u);

// Cosine-weighted hemisphere PDF. The closer the point/direction is to the top of the
// hemisphere, the higher the probability of sampling it will be. Theta is measured with
// respect to the axis of the hemisphere.
inline Float CosineHemispherePdf(Float cosTheta) {
    return cosTheta * InvPi;
}

// Samples the exponential distribution with parameter lambda using the Inverse Transform Method.
// Given a sample u ~ Unif(0,1), the random variable T = F^-1(u) has an exponential distribution,
// where the inverse of the exponential CDF
//
// F^-1(p) = -ln(1-p)/lambda
//
// is the quantile function of the exponential distribution, which gives the value of the exponential
// random variable such that the probability of the variable being less than or equal to that value
// equals the given probability p.
inline Float SampleExponential(Float u, Float lambda) {
    return -std::log(1 - u) / lambda;
}

// A piecewise-constant 1D distribution.
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

// A piecewise-constant 2D distribution.
class Distribution2D {
private:
    // Conditional density function p(u|v) of u given v for every value of v.
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;

    // Marginal density function p(v) for every value of v.
    std::unique_ptr<Distribution1D> pMarginal;

public:
    // The given array is a discretized 2D function, f(u,v).
    Distribution2D(const Float *f, int nu, int nv) {
        for (int v = 0; v < nv; ++v) {
            // Compute conditional sampling distribution for v. Note that the conditional
            // density function of u given v is based on the vth row.
            pConditionalV.emplace_back(new Distribution1D(&f[v * nu], nu));
        }

        // Compute marginal sampling distribution for v.
        std::vector<Float> marginalFunctions;
        for (int v = 0; v < nv; ++v) {
            marginalFunctions.push_back(pConditionalV[v]->definiteIntegral);
        }
        pMarginal.reset(new Distribution1D(&marginalFunctions[0], nv));
    }

    // Samples the 2D (u,v) distribution using a 2D uniformly distributed random variable, u.
    Point2f SampleContinuous(const Point2f &u, Float *pdf) const {
        Float pdfs[2];
        int v;

        // The 2nd component of the uniform bivariate sample is used to obtain a value of v
        // uniformly at random, as well as its marginal density p(v).
        Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);

        // The 1st component of the uniform bivariate sample is used to sample the conditional
        // distribution of u given v. 
        Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);

        *pdf = pdfs[0] * pdfs[1];
        return Point2f(d0, d1);
    }

    // PDF of the distribution. Gives the probability of sampling p (as sampled with SampleContinuous).
    Float Pdf(const Point2f &p) const {
        int uIndex = Clamp(int(p[0] * pConditionalV[0]->Count()), 0, pConditionalV[0]->Count() - 1);
        int vIndex = Clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[vIndex]->f[uIndex] / pMarginal->definiteIntegral;
    }
};

// Sample weighting function for multiple importance sampling, used to reduce variance
// of the Monte Carlo estimator for the integral of a product of functions f and g caused
// by mismatches between the integrand and the PDFs in the estimates f(X)g(X)/fPdf(X) or 
// f(X)g(X)/gPdf(X).
inline Float BalanceHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}

// Sample weighting function for multiple importance sampling, used to reduce variance
// of the Monte Carlo estimator for the integral of a product of functions f and g caused
// by mismatches between the integrand and the PDFs in the estimates f(X)g(X)/fPdf(X) or 
// f(X)g(X)/gPdf(X).
// 
// Reduces variance even more than the balance heuristic. Here, the parameter beta is
// always 2.
inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    Float f = nf * fPdf;
    Float g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

#endif // CPBRT_CORE_SAMPLING_H