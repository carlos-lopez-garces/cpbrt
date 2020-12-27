#include <memory>

#include "geometry.h"
#include "spectrum.h"

// TODO: implement.
std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name, Point2i *resolution);

// TODO: implement.
void WriteImage(
    const std::string &name,
    const Float *rgb,
    const Bounds2i &outputBounds,
    const Point2i &totalResolution
);