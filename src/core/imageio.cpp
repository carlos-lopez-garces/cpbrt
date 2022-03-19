#include "cpbrt.h"
#include "ext/lodepng.h"
#include "fileutil.h"
#include "imageio.h"
#include "spectrum.h"

// #include <ImfRgba.h>
// #include <ImfRgbaFile.h>

#define BUFFER_SIZE 80

static CPBRT_CONSTEXPR bool hostLittleEndian =
#if defined(__BYTE_ORDER__)
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    true
  #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    false
  #else
    #error "__BYTE_ORDER__ defined but has unexpected value"
  #endif
#else
  #if defined(__LITTLE_ENDIAN__) || defined(__i386__) || defined(__x86_64__) || \
      defined(_WIN32) || defined(WIN32)
    true
  #elif defined(__BIG_ENDIAN__)
    false
  #elif defined(__sparc) || defined(__sparc__)
    false
  #else
    #error "Can't detect machine endian-ness at compile-time."
  #endif
#endif
;

static inline int isWhitespace(char c) {
    return c == ' ' || c == '\n' || c == '\t';
}

static int readWord(FILE *fp, char *buffer, int bufferLength) {
    int n;
    int c;

    if (bufferLength < 1) return -1;

    n = 0;
    c = fgetc(fp);
    while (c != EOF && !isWhitespace(c) && n < bufferLength) {
        buffer[n] = c;
        ++n;
        c = fgetc(fp);
    }

    if (n < bufferLength) {
        buffer[n] = '\0';
        return n;
    }

    return -1;
}

static RGBSpectrum *ReadImagePFM(
    const std::string &filename, int *xres, int *yres
) {
    float *data = nullptr;
    RGBSpectrum *rgb = nullptr;
    char buffer[BUFFER_SIZE];
    unsigned int nFloats;
    int nChannels, width, height;
    float scale;
    bool fileLittleEndian;

    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) goto fail;

    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;

    if (strcmp(buffer, "Pf") == 0)
        nChannels = 1;
    else if (strcmp(buffer, "PF") == 0)
        nChannels = 3;
    else
        goto fail;

    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    width = atoi(buffer);
    *xres = width;

    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    height = atoi(buffer);
    *yres = height;

    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    sscanf(buffer, "%f", &scale);

    nFloats = nChannels * width * height;
    data = new float[nFloats];
    for (int y = height - 1; y >= 0; --y) {
        if (fread(&data[y * nChannels * width], sizeof(float),
                  nChannels * width, fp) != nChannels * width)
            goto fail;
    }

    fileLittleEndian = (scale < 0.f);
    if (hostLittleEndian ^ fileLittleEndian) {
        uint8_t bytes[4];
        for (unsigned int i = 0; i < nFloats; ++i) {
            memcpy(bytes, &data[i], 4);
            std::swap(bytes[0], bytes[3]);
            std::swap(bytes[1], bytes[2]);
            memcpy(&data[i], bytes, 4);
        }
    }
    if (std::abs(scale) != 1.f)
        for (unsigned int i = 0; i < nFloats; ++i) data[i] *= std::abs(scale);

    rgb = new RGBSpectrum[width * height];
    if (nChannels == 1) {
        for (int i = 0; i < width * height; ++i) rgb[i] = RGBSpectrum(data[i]);
    } else {
        for (int i = 0; i < width * height; ++i) {
            Float frgb[3] = {data[3 * i], data[3 * i + 1], data[3 * i + 2]};
            rgb[i] = RGBSpectrum::FromRGB(frgb);
        }
    }

    delete[] data;
    fclose(fp);
    return rgb;

fail:
    Error("Error reading PFM file \"%s\"", filename.c_str());
    if (fp) fclose(fp);
    delete[] data;
    delete[] rgb;
    return nullptr;
}

std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name, Point2i *resolution) {
    if (HasExtension(name, ".pfm")) {
        return std::unique_ptr<RGBSpectrum[]>(
            ReadImagePFM(name, &resolution->x, &resolution->y)
        );
    }

    Error("Unable to load image stored in format \"%s\" for filename \"%s\".",
          strrchr(name.c_str(), '.') ? (strrchr(name.c_str(), '.') + 1)
                                     : "(unknown)",
          name.c_str());
    return nullptr;
}

// static void WriteImageEXR(const std::string &name, const Float *pixels,
//                           int xRes, int yRes, int totalXRes, int totalYRes,
//                           int xOffset, int yOffset) {
//     using namespace Imf;
//     using namespace Imath;

//     Rgba *hrgba = new Rgba[xRes * yRes];
//     for (int i = 0; i < xRes * yRes; ++i)
//         hrgba[i] = Rgba(pixels[3 * i], pixels[3 * i + 1], pixels[3 * i + 2]);

//     // OpenEXR uses inclusive pixel bounds.
//     Box2i displayWindow(V2i(0, 0), V2i(totalXRes - 1, totalYRes - 1));
//     Box2i dataWindow(V2i(xOffset, yOffset),
//                      V2i(xOffset + xRes - 1, yOffset + yRes - 1));

//     try {
//         RgbaOutputFile file(name.c_str(), displayWindow, dataWindow,
//                             WRITE_RGB);
//         file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
//         file.writePixels(yRes);
//     } catch (const std::exception &exc) {
//         Error("Error writing \"%s\": %s", name.c_str(), exc.what());
//     }

//     delete[] hrgba;
// }

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
    } else if (HasExtension(name, ".exr")) {
        // WriteImageEXR(
        //     name, rgb, resolution.x, resolution.y, totalResolution.x, totalResolution.y, outputBounds.pMin.x, outputBounds.pMin.y
        // );
    } else {
        Error("Can't determine image file type from suffix of filename \"%s\"",
              name.c_str());
    }
}