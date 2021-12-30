#include "homogeneous.h"

Spectrum HomogeneousMedium::Tr(const Ray &ray, Sampler &sampler) const {
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}