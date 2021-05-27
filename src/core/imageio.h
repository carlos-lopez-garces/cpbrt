#ifndef CPBRT_CORE_IMAGEIO_H
#define CPBRT_CORE_IMAGEIO_H

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

#endif // CPBRT_CORE_IMAGEIO_H