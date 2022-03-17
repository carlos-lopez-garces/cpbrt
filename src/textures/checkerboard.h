#ifndef CPBRT_TEXTURES_CHECKERBOARD_H
#define CPBRT_TEXTURES_CHECKERBOARD_H

#include "core/cpbrt.h"
#include "core/texture.h"
#include "core/interaction.h"
#include "core/paramset.h"
#include "core/spectrum.h"

enum class AAMethod { 
    // Simple point sampling. No antialiasing.
    None,

    // Closed-form box filter.
    ClosedForm
};

template <typename T> class Checkerboard2DTexture : public Texture<T> {
private:
    // Maps the intersection point in the SurfaceInteraction to a 2D (s,t) texture
    // coordinate.
    std::unique_ptr<TextureMapping2D> mapping;

    // The checkerboard alternates between these 2 textures.
    const std::shared_ptr<Texture<T>> tex1, tex2;

    // Antialiasing method.
    const AAMethod aaMethod;

public:
    Checkerboard2DTexture(
        std::unique_ptr<TextureMapping2D> mapping,
        const std::shared_ptr<Texture<T>> &tex1,
        const std::shared_ptr<Texture<T>> &tex2,
        AAMethod aaMethod
    ) : mappping(mapping), tex1(tex1), tex2(tex2), aaMethod(aaMethod) {}

    T Evaluate(const SurfaceInteraction &si) const {
        // Partial-s-partial-x and partial-t-partial-x.
        Vector2f dstdx;
        // Partial-s-partial-y and partial-t-partial-y.
        Vector2f dstdy;

        // Maps the intersection point (which is in world space) to a texture
        // (s,t) coordinate.
        Point2f st = mapping->Map(si, &dstdx, &dstdy);

        if (aaMethod == AAMethod::None) {
            // Point sampling (no antialiasing): if floor(s) + floor(t) is even, sample
            // the first texture; if it's odd, sample the second texture.
            if (((int) std::floor(st[0]) + (int) std::floor(st[1])) % 2 == 0) {
                return tex1->Evaluate(si);
            }
            return tex2->Evaluate(si);
        } else {
            // Compute closed-form box-filtered value.

            // Define an axis-aligned bounding box around (s,t) of width 2ds and
            // height 2dt.
            Float ds = std::max(std::abs(dstdx[0]), std::abs(dstdy[0]));
            Float dt = std::max(std::abs(dstdx[1]), std::abs(dstdy[1]));
            Float s0 = st[0] - ds;
            Float s1 = st[0] + ds;
            Float t0 = st[1] - dt;
            Float t1 = st[1] + dt;

            // Each check of the checkerboard has unit side length. If the discretized
            // endpoints of the bounding box width and the discretized endpoints of its
            // height coincide, then the bounding box is contained inside a single check;
            // point-sample the check in that case.
            if (std::floor(s0) == std::floor(s1) 
                && std::floor(t0) == std::floor(t1)) {

                // Point sampling (no antialiasing): if floor(s) + floor(t) is even, sample
                // the first texture; if it's odd, sample the second texture. 
                if (((int) std::floor(st[0]) + (int) std::floor(st[1])) % 2 == 0) {
                    return tex1->Evaluate(si);
                }
                return tex2->Evaluate(si);
            }

            // The filter extent spans checks of the 2 types. Obtain the fraction of the
            // filters extent that covers each check type and use them to interpolate between
            // a sample from one subtexture and a sample from the other.
            //
            // Let c(x) be a 1D step function that maps x to 0 if x lies in check of type 1 and
            // to 1 if it lies in check of type 2, where -ds < x < ds (-dt < x < dt). Then integrate
            // c(x) in the interval (0, x) to compute the average of c(x), which gives the fraction of
            // the dominant check type (the check type that the box covers the most).
            //
            // Consider that the box filter may span multiple checks of both types.

            // Definite integral of 1D step function c(x) with (0, x) interval of integration.
            // Gives the average value of c(x).
            auto stepCofXIntegral = [](Float x) {
                return (int) std::floor(x/2)
                    + 2*std::max(x/2 - (int) std::floor(x/2) - (Float) 0.5, (Float) 0);
            };

            Float sAvg = (stepCofXIntegral(s1) - stepCofXIntegral(s0)) / (2*ds);
            Float tAvg = (stepCofXIntegral(t1) - stepCofXIntegral(t0)) / (2*dt);
            Float area2 = sAvg + tAvg - 2*sAvg*tAvg;
            if (ds > 1 || dt > 1) {
                area2 = .5f;
            }

            // Interpolate between the 2 subtextures using the fraction of the box filter
            // that covers each.
            return (1 - area2)*tex1->Evaluate(si) + area2*tex2->Evaluate(si);
        }
    }
};

Texture<Spectrum> *CreateCheckerboardSpectrumTexture(const Transform &tex2world, const TextureParams &tp);

#endif // CPBRT_TEXTURES_CHECKERBOARD_H