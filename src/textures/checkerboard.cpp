#include "textures/checkerboard.h"

Texture<Spectrum> *CreateCheckerboardSpectrumTexture(
    const Transform &tex2world, const TextureParams &tp
) {
    int dim = tp.FindInt("dimension", 2);
    if (dim != 2) {
        // TODO: add Checkerboard3DTexture.
        Error("%d dimensional checkerboard texture not supported", dim);
        return nullptr;
    }

    std::shared_ptr<Texture<Spectrum>> tex1 = tp.GetSpectrumTexture("tex1", 1.f);
    std::shared_ptr<Texture<Spectrum>> tex2 = tp.GetSpectrumTexture("tex2", 0.f);

    if (dim == 2) {
        // Mapping.
        std::unique_ptr<TextureMapping2D> map;
        std::string type = tp.FindString("mapping", "uv");
        if (type == "uv") {
            Float su = tp.FindFloat("uscale", 1.);
            Float sv = tp.FindFloat("vscale", 1.);
            Float du = tp.FindFloat("udelta", 0.);
            Float dv = tp.FindFloat("vdelta", 0.);
            map.reset(new UVMapping2D(su, sv, du, dv));
        } else {
            // TODO: add SphericalMapping2D, CylindricalMapping2D, and PlanarMapping2D.
            Error("2D texture mapping \"%s\" unknown", type.c_str());
            map.reset(new UVMapping2D);
        }

        // AAMethod.
        std::string aa = tp.FindString("aamode", "closedform");
        AAMethod aaMethod;
        if (aa == "none") {
            aaMethod = AAMethod::None;
        }
        else if (aa == "closedform") {
            aaMethod = AAMethod::ClosedForm;
        } else {
            Warning(
                "Antialiasing mode \"%s\" not understood by "
                "Checkerboard2DTexture; using \"closedform\"",
                aa.c_str());
            aaMethod = AAMethod::ClosedForm;
        }

        return new Checkerboard2DTexture<Spectrum>(std::move(map), tex1, tex2, aaMethod);
    } else {
        // TODO: add Checkerboard3DTexture.
      return nullptr;
    }
}