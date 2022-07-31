#ifndef CPBRT_CORE_INTERPOLATION_H
#define CPBRT_CORE_INTERPOLATION_H

// Computes the weights needed to obtain a spline-based interpolator of a function f and
// evaluate it at x, using 4 controls points. The interpolator is a polynomial of the form:
//
// p(x) = w0*f(x_-1) + w1*f(x_0) + w2*f(x_1) + w3*f(x_2).
//
// The offset is the index of the node that corresponds to the start endpoint x_i of the
// subinterval [x_i, x_i+1] that contains x. That's in 1D. In general, the nodes or points
// of the spline act as subinterval endpoints and the line segments between them as the
// subintervals.
bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset, Float *weights);

#endif // CPBRT_CORE_INTERPOLATION_H