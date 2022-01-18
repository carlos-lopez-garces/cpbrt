#include "lightdistribution.h"

UniformLightDistribution::UniformLightDistribution(const Scene &scene) {
    std::vector<Float> prob(scene.lights.size(), Float(1));
    distrib.reset(new Distribution1D(&prob[0], int(prob.size())));
}

const Distribution1D *UniformLightDistribution::Lookup(const Point3f &p) const {
    return distrib.get();
}

std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
    const std::string &name, const Scene &scene) {
    if (name == "uniform" || scene.lights.size() == 1) {
        return std::unique_ptr<LightDistribution>{
            new UniformLightDistribution(scene)
        };
    } else {
        Error(
            "Light sample distribution type \"%s\" unknown. Using \"uniform\".",
            name.c_str()
        );
        return std::unique_ptr<LightDistribution>{
            new UniformLightDistribution(scene)
        };
    }
}