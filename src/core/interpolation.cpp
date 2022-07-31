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
    *offset - nodeIdx - 1;
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