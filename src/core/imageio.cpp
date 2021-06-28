#include "cpbrt.h"
#include "ext/lodepng.h"
#include "imageio.h"
#include "fileutil.h"

void WriteImage(
    const std::string &name,
    const Float *rgb,
    const Bounds2i &outputBounds,
    const Point2i &totalResolution
) {
    Vector2i resolution = outputBounds.Diagonal();
    // TODO: add other types.
    if (HasExtension(name, ".png")) {
        // 8-bit formats; apply gamma
        Vector2i resolution = outputBounds.Diagonal();
        std::unique_ptr<uint8_t[]> rgb8(
            new uint8_t[3 * resolution.x * resolution.y]);
        uint8_t *dst = rgb8.get();
        for (int y = 0; y < resolution.y; ++y) {
            for (int x = 0; x < resolution.x; ++x) {
#define TO_BYTE(v) (uint8_t) Clamp(255.f * GammaCorrect(v) + 0.5f, 0.f, 255.f)
                dst[0] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 0]);
                dst[1] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 1]);
                dst[2] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 2]);
#undef TO_BYTE
                dst += 3;
            }
        }

        unsigned int error = lodepng_encode24_file(
            name.c_str(), rgb8.get(), resolution.x, resolution.y);
        if (error != 0)
            Error("Error writing PNG \"%s\": %s", name.c_str(),
                    lodepng_error_text(error));
    } else {
        Error("Can't determine image file type from suffix of filename \"%s\"",
              name.c_str());
    }
}