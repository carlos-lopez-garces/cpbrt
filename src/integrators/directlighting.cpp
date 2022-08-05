#include "directlighting.h"
#include "core/interaction.h"
#include "core/paramset.h"
#include "core/camera.h"
#include "core/film.h"

void DirectLightingIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    if (strategy == LightStrategy::UniformSampleAll) {
        // Compute number of samples to use for each light source.
        for (const auto &light : scene.lights) {
            nLightSamples.push_back(sampler.RoundCount(light->nSamples));
        }

        // Request samples for sampling all lights. Taking one sample from a light source
        // requires 2 2D [0,1) samples, one for sampling the surface of the light (see 
        // Shape::Sample) and one for sampling the BSDF at the point being shaded (see
        // BSDF::Sample_f). Taking many such samples (nLightSamples) from all the light 
        // sources requires 2 array-samples of nLightSamples 2D [0,1) samples each: the first
        // array-sample is for sampling the light sources' surfaces and the second one for the
        // BSDFs.
        //
        // Note that samples are requested for every light at every tracing depth level, i.e.
        // at every possible ray intersection point.
        for (int i = 0; i < maxDepth; ++i) {
            for (size_t j = 0; j < scene.lights.size(); ++j) {
                // For sampling points from the surface of light sources. 
                sampler.Request2DArray(nLightSamples[j]);
                // For sampling BSDFs at shaded points.
                sampler.Request2DArray(nLightSamples[j]);
            }
        }
    }
}

Spectrum DirectLightingIntegrator::Li(
    const RayDifferential &ray,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    int depth
) const {
    // Outgoing? radiance to compute.
    Spectrum L(0.f);

    SurfaceInteraction si;
    if (!scene.Intersect(ray, &si)) {
        // The ray escapes the scene bounds without having hit anything. Sample environment
        // lights, if there are any.
        for (const auto &light : scene.lights) {
            L += light->Le(ray);
        }

        // Tracing stops here.
        return L;
    }

    // Create and compute BSDF at the surface-ray intersection point and associate it with the
    // Interaction.
    si.ComputeScatteringFunctions(ray, arena);
    if (!si.bsdf) {
        // TODO: No BSDF, so the ray passes right through?
        return Li(si.SpawnRay(ray.d), scene, sampler, arena, depth);
    }

    // Outgoing direction as computed by the BSDF.
    Vector3f wo = si.wo;

    // If the intersected object is emissive (an area light, for example), it contributes to the
    // radiance carried by the ray.
    L += si.Le(wo);

    if (scene.lights.size() > 0) {
        // Estimate the value of the outgoing? radiance integral at the point of intersection.
        if (strategy == LightStrategy::UniformSampleAll) {
            // Take a number of samples from all of the light sources.
            L += UniformSampleAllLights(si, scene, arena, sampler, nLightSamples);
        } else {
            // Take a single sample from a single light source chosen uniformly at random.
            L += UniformSampleOneLight(si, scene, arena, sampler);
        }
    }

    if (depth+1 < maxDepth) {
        // Trace rays recursively for specular reflection and transmission. In general, the direct
        // lighting integrator estimates incident radiance using samples from light sources that
        // illuminate the surface directly. But in order for a surface with a (perfect) specular
        // BDRF to reflect the image of nearby objects, the integrator needs to trace and accumulate
        // radiance that bounces off of those objects in the direction of perfect specular reflection
        // toward the intersection point.

        // SpecularReflect calls Li recursively until depth reaches maxDepth.
        L += SpecularReflect(ray, si, scene, sampler, arena, depth);
        // TODO: implement.
        // L += SpecularTransmit(ray, si, scene, sampler, arena, depth);
    }

    return L;
}

DirectLightingIntegrator *CreateDirectLightingIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    LightStrategy strategy;
    std::string st = params.FindOneString("strategy", "all");
    if (st == "one")
        strategy = LightStrategy::UniformSampleOne;
    else if (st == "all")
        strategy = LightStrategy::UniformSampleAll;
    else {
        Warning(
            "Strategy \"%s\" for direct lighting unknown. "
            "Using \"all\".",
            st.c_str());
        strategy = LightStrategy::UniformSampleAll;
    }
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
    return new DirectLightingIntegrator(strategy, maxDepth, camera, sampler);
}