#include "texture.h"
#include "shape.h"

Point2f UVMapping2D::Map(const SurfaceInteraction &si, Vector2f *dstdx, Vector2f *dstdy) const {
    // Let T be the mapping (u,v) -> (s,t) as defined by the returned value of this Map function.
    // Its derivative is then:
    //
    // T'(u,v) = ((ds/du, ds/dv), (dt/du, dt/dv)) = ((su, 0), (0, sv))
    //
    // where su and sv are the scaling factors and are not to be confused with the s coordinate;
    // likewise, du and dv in the returned value of this Map function are the shifting deltas and
    // are not to be confused with the du and dv partial differentials used above.
    //
    // In order to match the texture function sampling rate to the image sampling rate, we need
    // to make u and v functions of x and y pixel coordinates: u(x,y) and v(x,y).
    //
    // T thus becomes a vector-valued function of the vector (x,y): T: (u(x,y),v(x,y)) -> (s,t).
    // Its derivative is then:
    //
    // T'(u(x,y), v(x,y)) = (dst/dx, dst/dy)
    //
    // where the dstdx and dstdy partial derivatives are given by the chain rule:
    //
    // dst/dx = (ds/dx, dt/dx) = (du/dx*ds/du + dv/dx*ds/dv, du/dx*dt/du + dv/dx*dt/dv)
    //                         = (du/dx*su    + dv/dx*0,     du/dx*0     + dv/dx*sv)
    //                         = (du/dx*su, dv/dx*sv)
    //
    // dst/dy = (ds/dy, dt/dy) = (du/dy*ds/du + dv/dy*ds/dv, du/dy*dt/du + dv/dy*dt/dv)
    //                         = (du/dy*su    + dv/dy*0,     du/dy*0     + dv/dy*sv)
    //                         = (du/dy*su, dv/dy*sv)
    *dstdx = Vector2f(si.dudx * su, si.dvdx * sv);
    *dstdy = Vector2f(si.dudy * su, si.dvdy * sv);

    // Scale and shift the (u,v) coordinate.
    return Point2f(su * si.uv[0] + du, sv * si.uv[1] + dv);
}