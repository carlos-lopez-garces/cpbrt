#ifndef CPBRT_CORE_TEXTURE_H
#define CPBRT_CORE_TEXTURE_H

#include "cpbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

class TextureMapping2D {
public:
    // Maps the intersection point in the SurfaceInteraction to a 2D (s,t) texture
    // coordinate. dst/dx and dst/dy are the partial derivatives of this mapping
    // with respect to a pixel-coordinate change in raster space; they tell the
    // texture function sampling rate that corresponds to the image function sampling
    // rate carried by the SurfaceInteraction.
    virtual Point2f Map(
        const SurfaceInteraction &si,
        Vector2f *dstdx,
        Vector2f *dstdy
    ) const = 0;
};

class UVMapping2D : public TextureMapping2D {
private:
    // Scaling factors for (u,v) coordinates applied when mapping a point on a surface to
    // a texture (s,t)=(u,v) coordinate.
    const Float su;
    const Float sv;

    // Shifting deltas for (u,v) coordinates applied when mapping a point on a surface to
    // a texture (s,t)=(u,v) coordinate.
    const Float du;
    const Float dv;

public:
    UVMapping2D(Float su = 1, Float sv = 1, Float du = 0, Float dv = 0);

    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx, Vector2f *dstdy) const; 
};

template <typename T> class Texture {
public:
    // Evaluates the texture at the point mapped to the intersection point.
    virtual T Evaluate(const SurfaceInteraction &si) const = 0;
};

// The Lanczos windowed sinc function limits the extent of the sinc function to a
// number of cycles controlled by the tau parameter.
Float LanczosFilter(Float coordinate, Float tau = 2);

#endif // CPBRT_CORE_TEXTURE_H