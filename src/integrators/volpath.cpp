#include "volpath.h"
#include "core/error.h"
#include "core/paramset.h"

void VolPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
}

VolPathIntegrator *CreateVolPathIntegrator(
    const ParamSet &params, 
    std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera
) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy = params.FindOneString("lightsamplestrategy", "spatial");
    return new VolPathIntegrator(maxDepth, camera, sampler, pixelBounds, rrThreshold, lightStrategy);
}