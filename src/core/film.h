#include <memory>
#include <string>

#include "filter.h"
#include "geometry.h"
#include "parallel.h"

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
    };

    // All the pixels of the image, possibly cropped.
    std::unique_ptr<Pixel[]> pixels;

    // Stores filter weights computed at construction-time for 16x16 uniformly distributed
    // sample points over the 2D extent of the filter.
    static constexpr int filterTableWidth = 16;
    Float filterTable[filterTableWidth * filterTableWidth];

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

    // ?
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

    // Computes the bounds of the rectangular sampling region for the (possibly cropped) image.
    // These bounds extend beyond the film or crop window area to accomodate the entire extent
    // of the filters centered at the boundary pixels or the pixels near them. 
    Bounds2i GetSampleBounds() const;

    // Computes the physical bounds of the film.
    Bounds2f GetPhysicalExtent() const;
};