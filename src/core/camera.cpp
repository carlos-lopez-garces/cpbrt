#include "camera.h"

Float Camera::GenerateRayDifferential(
    const CameraSample &sample,
    RayDifferential *rd
) const {
    Float wt = GenerateRay(sample, rd);

    // Compute ray differential using sample shifted one pixel in the X direction on
    // the film plane; the point on the lens and the origin of the ray may change as
    // a result according to the Camera implementation.
    CameraSample xShiftedSample = sample;
    xShiftedSample.pFilm.x++;
    Ray rx;
    Float wtx = GenerateRay(xShiftedSample, &rx);
    if (wtx == 0) {
        return 0;
    }
    rd->rxOrigin = rx.o;
    rd->rxDirection = rx.d;

    // Compute ray differential using sample shifted one pixel in the Y direction on
    // the film plane.
    CameraSample yShiftedSample = sample;
    yShiftedSample.pFilm.y++;
    Ray ry;
    Float wty = GenerateRay(yShiftedSample, &ry);
    if (wty == 0) {
        return 0;
    }
    rd->ryOrigin = ry.o;
    rd->ryDirection = ry.d;

    rd->hasDifferentials = true;
    return wt;
};