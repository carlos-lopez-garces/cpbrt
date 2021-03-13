#include <memory>

#include "scene.h"

class Integrator {
public:
    virtual void Render(const Scene &scene) = 0;
};

class SamplerIntegrator : public Integrator {
private:
    // Supplies sample points on the image at which to compute radiance.
    std::shared_ptr<Sampler> sampler;

protected:
    std::shared_ptr<const Camera> camera;

public:
    SamplerIntegrator(
        std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler
    ) : camera(camera), sampler(sampler) {}

    virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
};