#include "film.h"
#include "geometry.h"
#include "integrator.h"
#include "memory.h"

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