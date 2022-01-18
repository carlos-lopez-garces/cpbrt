#include "volpath.h"

void VolPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
}