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

    void Render(const Scene &scene);

    virtual void Preprocess(const Scene &scene, Sampler &sampler) {}

    // Samples the incident radiance function using the given ray (and differentials).
    virtual Spectrum Li(
        const RayDifferential &ray,
        const Scene &scene,
        Sampler &sampler,
        MemoryArena &arena,
        // Number of ray bounces from the camera that have occurred up until this call.
        int depth = 0
    ) const = 0;
};