#include "path.h"
#include "core/camera.h"
#include "core/film.h"
#include "core/interaction.h"
#include "core/paramset.h"
#include "core/scene.h"

Spectrum PathIntegrator::Li(
    const RayDifferential &r,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    int depth
) const {
    // Running total of the sum of path radiance contributions SUM(P(p)).
    Spectrum L(0.f);

    // Evaluation of the path's throughput function (of some of its factors, rather: when
    // the PDF over area in the light transport equation is expressed in terms of solid angle,
    // beta is the resulting expression of the throughput factor).
    Spectrum beta(1.f);

    RayDifferential ray(r);

    bool specularBounce = false;

    // Path sampling.
    for (int bounces = 0; ; ++bounces) {
        // Find next path vertex and accumulate contribution.

        // Intersect ray with scene to find next path vertex.
        SurfaceInteraction si;
        bool foundIntersection = scene.Intersect(ray, &si);

        // Possibly add emitted light at intersection.
        if (bounces == 0 || specularBounce) {
            // Add emitted light at path vertex or from the environment.
            //
            // When bounces = 0, this segment of the path starts directly at the camera.
            //
            // When specularBounce = true, the last segment of the path ended at a surface of
            // specular BSDF.

            if (foundIntersection) {
                // When bounces = 0, the length of the current path is i=2. In the general case,
                // radiance is only sampled at the i+1th vertex of the path (which is deliberately
                // a light source). But if the second vertex of the path (the first is on the camera)
                // is itself on a light source or emissive object, we can't ignore it, so in this
                // case we sample radiance at the ith vertex.
                L += beta * si.Le(-ray.d);
            } else {
                // The camera ray escaped out into the scene. Add the radiance contributions of
                // lights in the scene.
                for (const auto &light : scene.lights) {
                    L += beta * light->Le(ray);
                }
            }
        }

        if (!foundIntersection || bounces >= maxDepth) {
            // When no intersection was found, the ray escaped out into the environment.
            // Reaching the established maximum number of bounces also terminates path sampling.
            break;
        }

        // Compute BSDFs and skip over medium boundaries.
        si.ComputeScatteringFunctions(ray, arena, true);
        if (!si.bsdf) {
            // The boundaries between participating media with equal incides of refraction are
            // represented in the scene as surfaces with no BSDFs. Since their indices of refraction
            // are the same, neither the direction of the ray nor the radiance that it carries back
            // get altered. So just prolong the ray.
            ray = si.SpawnRay(ray.d);

            // These "intersections" don't count.
            bounces--;
            continue;
        }

        // Place the i+1th vertex of the path at a light source by sampling a point on one of them.
        // Compute the radiance contribution of the ith vertex (the current intersection) as a resut
        // of direct lighting from the chosen light source.
        L += beta * UniformSampleOneLight(si, scene, arena, sampler);

        // Computations for the path of length i finish here. The remainder of the iteration prepares
        // the next path, the path of length i+1. (Paths are not sampled from scratch; they are always
        // an extension of the previous one, so values like path throughput beta and total path radiance
        // L can be reused.)  

        // Sample the BSDF at the ith vertex to obtain a direction in which to extend the current path
        // of length i to obtain the next path of length i+i. Update path throughput beta.
        Vector3f wo = -ray.d;
        Vector3f wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = si.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) {
            break;
        }
        beta *= f * AbsDot(wi, si.shading.n) / pdf;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        ray = si.SpawnRay(wi);

        // TODO: Account for subsurface scattering, if applicable.

        if (bounces > 3) {
            // Terminate path probabilistically using Russian roulette.

            // The probability of termination q is inversely proportional to the current path throughput
            // beta. The smaller beta is at this vertex of the path, the smaller the contributions of
            // subsequent vertices are made (because light coming from subsequent vertices will necessarily
            // pass through this vertex and be subject to its BSDF; see the beta multiplication). A higher
            // probability of termination is desirable because path samples whose estimates contribute little
            // introduce variance in the total estimate of L.
            Float q = std::max((Float) 0.05, 1 - beta.y());
            if (sampler.Get1D() < q) {
                // Terminate.
                break;
            }

            // Weighting the path throughput beta compensates for the bias introduced by skipping path
            // samples when paths are terminated. 
            beta /= 1 - q;
        }
    }

    return L;
}

// TODO: explain.
PathIntegrator *CreatePathIntegrator(
    const ParamSet &params,
    std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera
) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new PathIntegrator(maxDepth, camera, sampler);
}