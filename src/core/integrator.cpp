#include "integrator.h"
#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "camera.h"

Spectrum UniformSampleAllLights(
    const Interaction &it,
    const Scene &scene,
    MemoryArena &arena,
    Sampler &sampler,
    const std::vector<int> &nLightSamples,
    bool handleMedia
) {
    Spectrum L(0.f);

    for (size_t j = 0; j < scene.lights.size(); ++j) {
        // Accumulate contribution of jth light to L.
        const std::shared_ptr<Light> &light = scene.lights[j];
        int nSamples = nLightSamples[j];

        // For sampling points from the surfaces of the light source.
        const Point2f *uLightArray = sampler.Get2DArray(nSamples);

        // For sampling the BSDF at the Interaction point for each of the directions
        // corresponding to each of the light source point samples.
        const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);

        // Evaluate the Monte Carlo estimator. By the law of large numbers, the expected
        // value of the estimator is the value of the direct lighting outgoing radiance integral.
        // The actual value is, of course, an approximation, because we only evaluate the
        // estimator using the values of n=nSamples uniform random variables.
        if (!uLightArray || !uScatteringArray) {
            // The preallocated array-samples have been exhausted, so we can't take the
            // number of samples that we wanted from all of the light sources. Take a single
            // point sample from a single light source instead.
            Point2f uLight = sampler.Get2D();
            Point2f uScattering = sampler.Get2D();

            // Evaluate the outgoing radiance / scattering equation for the incident direction
            // formed between the Interaction point and the random sample point on this light source.
            L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler, arena, handleMedia);
        } else {
            // Incident radiance coming from this light source.
            Spectrum Ld(0.f);

            // Evaluate and sum the outgoing radiance / scattering equation for the incident direction
            // formed between the Interaction point and each of the random sample points on this light
            // source.
            for (int k = 0; k < nSamples; ++k) {
                Ld += EstimateDirect(
                    it, uScatteringArray[k], *light, uLightArray[k], scene, sampler, arena, handleMedia
                );
            }

            // This is the value of the Monte Carlo estimator: the average outgoing radiance
            // corresponding to the n=nSamples incident directions as reflected by the BSDF.
            L += Ld / nSamples;
        }
    }

    return L;
}

Spectrum UniformSampleOneLight(
    const Interaction &it,
    const Scene &scene,
    MemoryArena &arena,
    Sampler &sampler,
    // Whether to account for the effects of volumetric attenuation.
    bool handleMedia
) {
    // Randomly choose single light to sample.
    int nLights = int(scene.lights.size());
    if (nLights == 0) {
        return Spectrum(0.f);
    }
    int lightNum = std::min((int) (sampler.Get1D() * nLights), nLights - 1);
    const std::shared_ptr<Light> &light = scene.lights[lightNum];

    // For sampling a point from the surface of the light source.
    Point2f uLight = sampler.Get2D();

    // For sampling the BSDF at the Interaction point for direction corresponding to light
    // source point sample.
    Point2f uScattering = sampler.Get2D();

    // Evaluate the outgoing radiance / scattering equation for this one direction. It can be
    // proven that the expected value of the sum of the scattering function evaluated for every
    // light source is equal to the value of the function evaluated just for one and multiplied
    // by the number of them.
    return (Float) nLights 
        * EstimateDirect(it, uScattering, *light, uLight, scene, sampler, arena, handleMedia);
}

Spectrum EstimateDirect(
    const Interaction &it,
    const Point2f &uScattering,
    const Light &light,
    const Point2f &uLight,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    bool handleMedia,
    bool specular
) {
    BxDFType bsdfFlags = specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);

    Spectrum Ld(0.f);

    // The scattering equation is a product of functions: the BRDF f and the direct radiance Ld from
    // the light. A distribution in the shape of this product is difficult to obtain. And if only the 
    // PDF of one of the 2 functions were to be used to sample the incident direction, the other function
    // would be sampled badly with it and the estimate would have high variance.
    //
    // Multiple importance sampling allows us to sample with the individual PDFs of the BRDF and of
    // Ld with low variance. 1 sample will be taken using the BRDF's PDF and 1 using the light source's
    // PDF, and their contributions to the direct radiance estimate will be weighted by a special 
    // weighting function. 

    // Sample light source with multiple importance sampling. The direction of incidence wi will be
    // determined from the light source's point distribution.
    Vector3f wi;
    Float lightPdf = 0;
    Float scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
    if (lightPdf > 0 && !Li.IsBlack()) {
        // Compute BSDF or phase function's value for light sample.
        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            // Evaluate BSDF for sampled incident direction.
            const SurfaceInteraction &si = (const SurfaceInteraction &) it;
            // Radiance is measured with respect to an area differential that is orthogonal
            // to the direction of incidence. To actually place this area differential on the
            // surface, the scattering equation includes the cosine of the theta angle as a factor,
            // measured from the surface normal to the direction of incidence). 
            f = si.bsdf->f(si.wo, wi, bsdfFlags) * AbsDot(wi, si.shading.n);
            scatteringPdf = si.bsdf->Pdf(si.wo, wi, bsdfFlags);
        } else {
            // TODO: Evaluate phase function for sampled incident direction. This is when the
            // Interaction is not a surface, but participating media.
        }

        if (!f.IsBlack()) {
            // The surface or medium reflects light.

            // Compute effect of visibility for light source sample.
            if (handleMedia) {
                Li *= visibility.Tr(scene, sampler);
            } else if (!visibility.Unoccluded(scene)) {
                // The light source doesn't illuminate the surface from the sampled direction.
                Li = Spectrum(0.f);
            }

            // Add light's contribution to reflected radiance.
            // TODO: This Li.IsBlack conditional seems redundant.
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags)) {
                    // The light source is described by a delta distribution, that is, it illuminates
                    // from a single direction with probability 1 and, therefore, the sample introduces
                    // no variance. Since there's no variance to reduce, there's no need for weighting
                    // the sample like multiple importance sampling does, so compute the standard Monte
                    // Carlo estimate for it and add its contribution.
                    Ld += f * Li / lightPdf;
                } else {
                    // Weight the sample and compute the MIS Monte Carlo estimate for it and add its
                    // contribution.
                    Float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                    Ld += f * Li * weight / lightPdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling. The direction of incidence wi will be determined
    // by the surface's (or medium's) BSDF.
    //
    // The surface's (or medium's) BSDF is not sampled when the light source is described by a delta
    // distribution, that is, when it illuminates from a single direction with probability 1: it would
    // be nearly impossible that a direction determined by sampling the BSDF correspond to this one
    // delta direction. (Radiance contributions from delta lights can only be obtained by sampling the
    // light source and has been done above already.) 
    if (!IsDeltaLight(light.flags)) {
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction()) {
            // Sample incident direction for surface interactions.
            BxDFType sampledType;
            const SurfaceInteraction &si = (const SurfaceInteraction &) it;
            f = si.bsdf->Sample_f(si.wo, &wi, uScattering, &scatteringPdf, bsdfFlags, &sampledType);
            // Radiance is measured with respect to an area differential that is orthogonal
            // to the direction of incidence. To actually place this area differential on the
            // surface, the scattering equation includes the cosine of the theta angle as a factor,
            // measured from the surface normal to the direction of incidence. 
            f *= AbsDot(wi, si.shading.n);
            sampledSpecular = sampledType & BSDF_SPECULAR;
        } else {
            // TODO: Sample scattered direction for medium interactions.
        }

        if (!f.IsBlack() && scatteringPdf > 0) {
            // Account for light contributions along sampled direction wi.
            //
            // If the BSDF is specular, wi was sampled with probability 1 (this incident direction
            // wi corresponds to the given scattered direction wo always), so its MIS weight must be
            // 1. The PDF already has the shape of the BSDF: delta?
            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);
                if (lightPdf == 0) {
                    // As computed with the light source sample only.
                    return Ld;
                }
                // Obtain the MIS weight of this sample. 
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }

            // Does the given light source illuminate this surface (medium) from the sampled incident
            // direction wi? Find intersection and compute transmittance.
            SurfaceInteraction lightIntersection;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction 
                = handleMedia ? scene.IntersectTr(ray, sampler, &lightIntersection, &Tr) 
                              : scene.Intersect(ray, &lightIntersection);

            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIntersection.primitive->GetAreaLight() == &light) {
                    // The given light source was intersected by a ray going in the sampled incident
                    // direction wi. Sample its emitted radiance.
                    Li = lightIntersection.Le(-wi);
                }
            } else {
                // The given light source doesn't have geometry in the scene (it is an infinite area
                // environment light, for example). Sample its emitted radiance.
                // TODO: What if it's occluded?
                Li = light.Le(ray);
            }

            // Compute the MIS Monte Carlo estimate for the BSDF (medium) sample and its contribution.
            if (!Li.IsBlack()) {
                Ld += f * Li * Tr * weight / scatteringPdf;
            }
        }
    }

    return Ld;
}

void SamplerIntegrator::Render(const Scene &scene) {
    Preprocess(scene, *sampler);

    // The bounds of the (possibly cropped) image with extra room to accomodate for the
    // extent of the reconstruction filters centered at pixels at or near the image boundaries.
    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();

    // The image will be broken up into tiles of 16x16 pixels.
    const int tileSize = 16;
    Point2i nTiles(
        (sampleExtent.x + tileSize - 1) / tileSize,
        (sampleExtent.y + tileSize - 1) / tileSize
    );

    // Render image tiles in parallel.
    ParallelFor2D(
        [&](Point2i tile) {
            // Render section of image corresponding to this tile.

            // Allocate MemoryArena for tile. A MemoryArena is a sizable chunk of memory
            // preallocated from the standard library memory manager and managed by a custom
            // allocator. Allocations are cheap because they are simple pointer increments, which
            // is great for making a large number of small allocations.
            MemoryArena arena;

            // Get sampler instance for tile. Samplers must not be shared across threads and it
            // is crucial that they use a unique seed so that they don't use the same sequences
            // of random numbers.
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile.
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            // If the image resolution (or crop window) is not a multiple of 16x16 tile size, cap
            // the far endpoint.
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            // If the image resolution (or crop window) is not a multiple of 16x16 tile size, cap
            // the far endpoint.
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

            // Get FilmTile for tile.
            std::unique_ptr<FilmTile> filmTile = camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them.
            for (Point2i pixel : tileBounds) {
                tileSampler->StartPixel(pixel);

                do {
                    // Initialize CameraSample for current sample.
                    CameraSample cameraSample = tileSampler->GetCameraSample(pixel);

                    // Generate camera ray for current sample.
                    RayDifferential ray;
                    Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(1 / std::sqrt(tileSampler->samplesPerPixel));

                    // Evaluate radiance along camera ray and arriving at the ray's origin.
                    Spectrum L(0.f);
                    if (rayWeight > 0) {
                        L = Li(ray, scene, *tileSampler, arena);
                    }

                    // Add camera ray's contribution to image.
                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free MemoryArena memory from computing image sample value.
                    arena.Reset();
                } while (tileSampler->StartNextSample());
            }

            // Merge image tile into Film.
            camera->film->MergeFilmTile(std::move(filmTile));
        },
        nTiles
    );

    // Save final image after rendering.
    camera->film->WriteImage();
}

Spectrum SamplerIntegrator::SpecularReflect(
    const RayDifferential &ray,
    const SurfaceInteraction &si,
    const Scene &scene,
    Sampler &sampler,
    MemoryArena &arena,
    int depth
) const {
    // Compute specular reflection direction wi and BSDF value. The outgoing direction wo
    // is the negated ray (camera or bounced/secondary) vector, and we want to find the incident
    // direction wi that reflects into wo.
    Vector3f wo = si.wo;
    Vector3f wi;
    Float pdf;
    // Only interested in evaluating specular BRDFs.
    BxDFType type = BxDFType(BSDF_REFLECTION |  BSDF_SPECULAR);
    Spectrum f = si->bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

    // Return contribution of specular reflection.
    const Normal3f &ns = si.shading.n;
    if (pdf > 0 && !f.IsBlack() && AbsDot(wi, ns) != 0) {
        // Compute ray differential for specular reflection. For antialiasing textures that get sampled
        // by a specular reflection ray.
        // TODO: explain.
        RayDifferential rd = si.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = si.p + si.dpdx;
            rd.ryOrigin = si.p + si.dpdy;

            // Compute differential reflected directions. 
            Normal3f dndx = si.shading.dndu * si.dudx + si.shading.dndv * si.dvdx;
            Normal3f dndy = si.shading.dndu * si.dudy + si.shading.dndv * si.dvdy;
            Vector3f dwodx = -ray.rxDirection - wo;
            Vector3f dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);
            rd.rxDirection = wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
            rd.ryDirection = wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
        }

        // Compute sum term of Monte Carlo estimator of scattering equation. The product of the
        // BRDF f and the incident radiance Li gives the fraction of incident light that will get
        // reflected. The AbsDot(wi, ns) = cos(wi, ns) factor places the area differential on the
        // surface (the area differential dA is originally perpendicular to the wi solid angle).
        return f * Li(rd, scene, sampler, arena, depth+1) * AbsDot(wi, ns) / pdf;
    } else {
        return Spectrum(0.f);
    }
}