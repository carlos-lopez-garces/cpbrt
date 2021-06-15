#ifndef CPBRT_CORE_FILM_H
#define CPBRT_CORE_FILM_H

#include <memory>
#include <mutex>
#include <string>

#include "cpbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "parallel.h"

// Stores 2 running sums, the numerator and denominator of the quotient that defines the
// pixel filtering/reconstruction equation:
//
// SUM_i(f(x-x_i, y-y_i)w(x_i, y_i)L(x_i, y_i))
// --------------------------------------------
// SUM_i(f(x-x_i, y-y_i)
//
// where f is the filter, w is the weight that the camera assigns to the sample, and L is
// radiance value of the sample.
struct FilmTilePixel {
    // Converting constructor of Spectrum allows for initialization from float.
    Spectrum contribSum = 0.f;
    Float filterWeightSum = 0.f;
};

class Film {
private:
    struct Pixel {
        // CIE XYZ color as reconstructed from pixel samples.
        Float xyz[3] = {0, 0, 0};

        Float filterWeightSum = 0;

        // ?
        AtomicFloat splatXYZ[3];

        // To avoid cache line straddling.
        Float pad;

        Pixel() { xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0; }
    };

    // All the pixels of the image, possibly cropped.
    // Layout: [(x0, y0), ..., (xn-1, y0), (x0, y1), ..., (xn-1, y1), ...]
    std::unique_ptr<Pixel[]> pixels;

    // Stores filter weights computed at construction-time for 16x16 uniformly distributed
    // sample points over the 2D extent of the filter.
    static constexpr int filterTableWidth = 16;
    Float filterTable[filterTableWidth * filterTableWidth];

    std::mutex mutex;

public:
    // In pixels.
    const Point2i fullResolution;

    // The length of the diagonal of the film's physical area, in meters.
    const Float diagonal;

    std::unique_ptr<Filter> filter;

    // Image file.
    const std::string filename;

    // Bounding box of subset of pixels to render as determined by a crop window. 
    Bounds2i croppedPixelBounds;

    // A user-supplied factor that scales the final value of the pixel (after sampling,
    // filtering, conversion to RGB, and splatting) for the user's purpose.
    const Float scale;

    Film(
        const Point2i &resolution,
        // Specifies a subset of pixels to render.
        const Bounds2f &cropWindow,
        std::unique_ptr<Filter> filt,
        // In millimeters.
        Float diagonal,
        const std::string &filename,
        Float scale
    ) : fullResolution(resolution),
        diagonal(diagonal * 0.001),
        filter(std::move(filt)),
        filename(filename),
        scale(scale) {
        
        // Compute film image bounds.
        croppedPixelBounds = Bounds2i(
            // Pixel located at the upper-left corner of the crop window (the y-axis
            // increases downward in NDC space, so the upper pixels have a low-value
            // y-coordinate).
            Point2i(
                std::ceil(fullResolution.x * cropWindow.pMin.x),
                std::ceil(fullResolution.y * cropWindow.pMin.y)
            ),

            // Pixel located at the lower-right corner of the crop window.
            Point2i(
                std::ceil(fullResolution.x * cropWindow.pMax.x),
                std::ceil(fullResolution.y * cropWindow.pMax.y)
            )
        );

        // Allocate film image storage.
        pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);

        // Precompute filter weight table. Filters are not evaluated for every single sample
        // point at rendering-time. Instead, the filter weights are precomputed and stored
        // in a lookup table for 16x16 uniformly distributed sample points over the 2D extent
        // of the filter. Then, at rendering-time, actual sample points look up their weight
        // at the entry of the precomputed sample point that is closest to them.
        //
        // The filter is assumed to be symmetric in both dimensions, so the filter weights of
        // positive coordinates are also the weights of the corresponding reflected coordinates. 
        int offset = 0;
        for (int y = 0; y < filterTableWidth; ++y) {
            for (int x = 0; x < filterTableWidth; ++x, ++offset) {
                Point2f p;
                // Sample point coordinates are confined to [0.5, radius.x - 0.5] and
                // [0.5, radius.y - 0.5] and equally spaced.
                p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
                p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
                filterTable[offset] = filter->Evaluate(p);
            }
        }
    }

    Pixel &GetPixel(const Point2i &p) {
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) + (p.y - croppedPixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    // Computes the bounds of the rectangular sampling region for the (possibly cropped) image.
    // These bounds extend beyond the film or crop window area to accomodate the entire extent
    // of the filters centered at the boundary pixels or the pixels near them. 
    Bounds2i GetSampleBounds() const;

    // Computes the physical bounds of the film.
    Bounds2f GetPhysicalExtent() const;

    // The film is divided into nonoverlapping tiles that rendering threads can sample and filter
    // independently.
    std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);

    void MergeFilmTile(std::unique_ptr<FilmTile> tile);

    // ?
    void SetImage(const Spectrum *image) const;

    // ?
    void AddSplat(const Point2f &p, const Spectrum &v);

    void WriteImage(Float splatScale = 1);
};

class FilmTile {
private:
    // Discrete pixel area of the film that the film tile covers.
    const Bounds2i pixelBounds;
    
    // Pixels enclosed by the film tile.
    // Layout: [(x0, y0), ..., (xn-1, y0), (x0, y1), ..., (xn-1, y1), ...]
    std::vector<FilmTilePixel> pixels;

    const Vector2f filterRadius;
    const Vector2f invFilterRadius;
    const Float *filterTable;

    // Filter table width (either number of rows or columns of square table).
    const int filterTableSize;

public:
    FilmTile(
        const Bounds2i &pixelBounds,
        const Vector2f &filterRadius,
        const Float *filterTable,
        int filterTableSize
    ) : pixelBounds(pixelBounds),
        filterRadius(filterRadius),
        invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
        filterTable(filterTable),
        filterTableSize(filterTableSize) {
        
        pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
    }

    FilmTilePixel &GetPixel(const Point2i &p) {
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset = (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }

    void AddSample(const Point2f &pFilm, const Spectrum &L, Float sampleWeight = 1.f) {
        // Compute discrete bounds that enclose pixels that the sample contributes to
        // according to the filter's extent. The filter is placed at the center of the
        // discrete pixel, not at the continuous pixel coordinate of the sample.

        // Continuous-to-discrete pixel coordinate mapping is: d = c - 0.5.
        Point2f pFilmDiscrete = pFilm - Vector2f(0.5f, 0.5f);
        Point2i p0 = (Point2i) Ceil(pFilmDiscrete - filterRadius);
        Point2i p1 = (Point2i) Floor(pFilmDiscrete + filterRadius) + Point2i(1, 1);
        // Shave off the pixels that are outside the image.
        p0 = Max(p0, pixelBounds.pMin);
        p1 = Min(p1, pixelBounds.pMax);

        // Precompute x and y filter table offsets (array indices). The offset is y*filterTableSize + x,
        // that's how the table was built by the Film constructor: for (y...width) { for (x...width) },
        // where width = filterTableSize.
        int *ifx = ALLOCA(int, p1.x - p0.x);
        for (int x = p0.x; x < p1.x; ++x) {
            // x is the x coordinate of a discrete pixel that is within the reach of the filter
            // centered at the discrete pixel pFilmDiscrete that contains the input sample point
            // pFilm. The distance between the 2 determines the wanted filter weight table x offset.
            Float fx = std::abs((x - pFilmDiscrete.x) * invFilterRadius.x * filterTableSize);
            ifx[x - p0.x] = std::min((int) std::floor(fx), filterTableSize - 1);
        }
        int *ify = ALLOCA(int, p1.y - p0.y);
        for (int y = p0.y; y < p1.y; ++y) {
            Float fy = std::abs((y - pFilmDiscrete.y) * invFilterRadius.y * filterTableSize);
            ify[y - p0.y] = std::min((int) std::floor(fy), filterTableSize - 1);
        }

        // Loop over discrete pixel coordinates within filter extent and add sample to
        // pixel arrays.
        for (int y = p0.y; y < p1.y; ++y) {
            for (int x = p0.x; x < p1.x; ++x) {
                // Evaluate filter weight for sample coordinate (look it up, really).
                //
                // This is how the filter would be evaluated if we weren't precomputing them:
                // filterWeight = filter->Evaluate(Point2i(x - pFilmDiscrete.x, y - pFilmDiscrete.y));
                int offset = ify[y - p0.y] * filterTableSize + ifx[x - p0.x];
                Float filterWeight = filterTable[offset];

                // Update pixel values with filtered sample contribution.
                FilmTilePixel &pixel = GetPixel(Point2i(x, y));
                pixel.contribSum += L * sampleWeight * filterWeight;
                pixel.filterWeightSum += filterWeight;
            }
        }
    }

    Bounds2i GetPixelBounds() const {
        return pixelBounds;
    }
};

#endif // CPBRT_CORE_FILM_H