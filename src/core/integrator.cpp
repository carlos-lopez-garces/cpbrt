#include "film.h"
#include "geometry.h"
#include "integrator.h"
#include "memory.h"

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
            Spectrumd Ld(0.f)

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