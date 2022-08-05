#include "cpbrt.h"
#include "interpolation.h"

// Computes the weights needed to obtain a spline-based interpolator of a function f and
// evaluate it at x, using 4 controls points. The interpolator is a polynomial of the form:
//
// p(x) = w0*f(x_-1) + w1*f(x_0) + w2*f(x_1) + w3*f(x_2).
//
// The offset is the index of the node that corresponds to the start endpoint x_i of the
// subinterval [x_i, x_i+1] that contains x. That's in 1D. In general, the nodes or points
// of the spline act as subinterval endpoints and the line segments between them as the
// subintervals.
bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset, Float *weights) {
    if (!(x >= nodes[0] && x <= nodes[size-1])) {
        // x is outside of the domain of the function.
        return false;
    }

    // Search for the interval [x0, x1] that contains x.
    int nodeIdx = FindInterval(size, [&](int i) {
        return nodes[i] <= x;
    });
    // Index of the start endpoint (spline node) of the subinterval (spline segment).
    *offset = nodeIdx - 1;
    Float x0 = nodes[nodeIdx];
    Float x1 = nodes[nodeIdx+1];


    // Compute the value of a parameter t relative to the unit interval [0,1] and that
    // corresponds to x relative to f's domain interval.
    Float t = (x - x0) / (x1 - x0);
    Float t2 = t * t;
    Float t3 = t2 * t;

    // Compute node weights w1, w2, w3, and w4.

    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;

    if (nodeIdx > 0) {
        Float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[nodeIdx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        Float w0 = t3 - 2 * t2 + t;
        weights[0] = 0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    if (nodeIdx + 2 < size) {
        Float w3 = (t3 - t2) * (x1 - x0) / (nodes[nodeIdx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        Float w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0;
    }

    return true;
}

Float SampleCatmullRom2D(
    int size1, int size2,
    const Float *nodes1,
    const Float *nodes2,
    const Float *values,
    const Float *cdf,
    Float alpha,
    Float u,
    Float *fval,
    Float *pdf
) {
    int offset;
    Float weights[4];
    if (!CatmullRomWeights(size1, nodes1, alpha, &offset, weights)) {
        return 0;
    }

    auto interpolate = [&](const Float *array, int idx) {
        Float value = 0;
        for (int i = 0; i < 4; ++i) {
            if (weights[i] != 0) {
                value += array[(offset + i) * size2 + idx] * weights[i];
            }
        }
        return value;
    };

    Float maximum = interpolate(cdf, size2 - 1);
    u *= maximum;
    int idx = FindInterval(size2, [&](int i) { 
        return interpolate(cdf, i) <= u; 
    });

    Float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
    Float x0 = nodes2[idx], x1 = nodes2[idx + 1];
    Float width = x1 - x0;
    Float d0, d1;

    u = (u - interpolate(cdf, idx)) / width;

    if (idx > 0) {
        d0 = width * (f1 - interpolate(values, idx - 1)) / (x1 - nodes2[idx - 1]);
    } else {
        d0 = f1 - f0;
    }
    if (idx + 2 < size2) {
        d1 = width * (interpolate(values, idx + 2) - f0) / (nodes2[idx + 2] - x0);
    } else {
        d1 = f1 - f0;
    }

    Float t;
    if (f0 != f1) {
        t = (f0 - std::sqrt(std::max((Float)0, f0 * f0 + 2 * u * (f1 - f0)))) / (f0 - f1);
    } else {
        t = u / f0;
    }

    Float a = 0, b = 1, Fhat, fhat;
    while (true) {
        if (!(t >= a && t <= b)) {
            t = 0.5f * (a + b);
        }

        Fhat = t * (f0 +
                    t * (.5f * d0 +
                         t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
                              t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 +
               t * (d0 +
                    t * (-2 * d0 - d1 + 3 * (f1 - f0) +
                         t * (d0 + d1 + 2 * (f0 - f1))));

        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) {
            break;
        }

        if (Fhat - u < 0) {
            a = t;
        }
        else {
            b = t;   
        }

        t -= (Fhat - u) / fhat;
    }

    if (fval) *fval = fhat;
    if (pdf) *pdf = fhat / maximum;
    return x0 + width * t;
}

Float IntegrateCatmullRom(int n, const Float *x, const Float *values, Float *cdf) {
    Float sum = 0;
    cdf[0] = 0;
    for (int i = 0; i < n - 1; ++i) {
        Float x0 = x[i], x1 = x[i + 1];
        Float f0 = values[i], f1 = values[i + 1];
        Float width = x1 - x0;

        Float d0, d1;
        if (i > 0) {
            d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
        } else {
            d0 = f1 - f0;
        }
        if (i + 2 < n) {
            d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
        } else {
            d1 = f1 - f0;
        }

        sum += ((d0 - d1) * (1.f / 12.f) + (f0 + f1) * .5f) * width;
        cdf[i + 1] = sum;
    }
    return sum;
}