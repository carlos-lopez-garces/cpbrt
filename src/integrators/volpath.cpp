#include "volpath.h"
#include "core/error.h"
#include "core/paramset.h"

void VolPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {}

Spectrum VolPathIntegrator::Li(
    const RayDifferential &r,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    int depth
) const {
    Spectrum L(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
    int bounces;
    Float etaScale = 1;

    for (bounces = 0;; ++bounces) {
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

        MediumInteraction mi;
        if (ray.medium) beta *= ray.medium->Sample(ray, sampler, arena, &mi);
        if (beta.IsBlack()) break;

        if (mi.IsValid()) {
            if (bounces >= maxDepth) break;

            L += beta * UniformSampleOneLight(mi, scene, arena, sampler, true);

            Vector3f wo = -ray.d, wi;
            mi.phase->Sample_p(wo, &wi, sampler.Get2D());
            ray = mi.SpawnRay(wi);
            specularBounce = false;
        } else {
            if (bounces == 0 || specularBounce) {
                if (foundIntersection) {
                    L += beta * isect.Le(-ray.d);
                } else {
                    for (const auto &light : scene.infiniteLights) {
                        L += beta * light->Le(ray);
                    }
                }
            }

            if (!foundIntersection || bounces >= maxDepth) break;

            isect.ComputeScatteringFunctions(ray, arena, true);
            if (!isect.bsdf) {
                ray = isect.SpawnRay(ray.d);
                bounces--;
                continue;
            }

            L += beta * UniformSampleOneLight(isect, scene, arena, sampler, true);

            Vector3f wo = -ray.d, wi;
            Float pdf;
            BxDFType flags;
            Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);
            if (f.IsBlack() || pdf == 0.f) break;
            beta *= f * AbsDot(wi, isect.shading.n) / pdf;
            Assert(std::isinf(beta.y()) == false);
            specularBounce = (flags & BSDF_SPECULAR) != 0;
            if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
                Float eta = isect.bsdf->eta;
                etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
            }
            ray = isect.SpawnRay(wi);

            if (isect.bssrdf && (flags & BSDF_TRANSMISSION)) {
                // TODO: Importance sample the BSSRDF.
            }
        }

        Spectrum rrBeta = beta * etaScale;
        if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
            Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q) break;
            beta /= 1 - q;
            Assert(std::isinf(beta.y()) == false);
        }
    }

    return L;
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
    std::string lightStrategy = params.FindOneString("lightsamplestrategy", "uniform");
    return new VolPathIntegrator(maxDepth, camera, sampler, pixelBounds, rrThreshold, lightStrategy);
}