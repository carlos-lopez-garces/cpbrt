#ifndef CPBRT_CORE_PARAMSET_H
#define CPBRT_CORE_PARAMSET_H

#include <memory>
#include <vector>

#include "cpbrt.h"
#include "geometry.h"

template <typename T> struct ParamSetItem {
    const std::string name;
    const std::unique_ptr<T[]> values;
    const int nValues;
    // Set to true when the item is read. If the item goes unread, it means that
    // it wasn't consumed/used, presumably because the name was misspelled in the 
    // input file.
    mutable bool lookedUp = false;

    template <typename T> ParamSetItem(const std::string &name, std::unique_ptr<T[]> value, int nValues)
      : name(name), values(new T[nValues]), nValues(nValues) {
          std::copy(v, v + nValues, values.get());
    }
};

class ParamSet {
public:

    // Return a supplied default value when the parameter is not found.
    bool FindOneBool(const std::string &name, bool defaultValue) const;
    int FindOneInt(const std::string &name, int defaultValue) const;
    Float FindOneFloat(const std::string &name, Float defaultValue) const;
    Point2f FindOnePoint2f(const std::string &name, const Point2f &defaultValue) const;
    Vector2f FindOneVector2f(const std::string &name, const Vector2f &defaultValue) const;
    Point3f FindOnePoint3f(const std::string &name, const Point3f &defaultValue) const;
    Vector3f FindOneVector3f(const std::string &name, const Vector3f &defaultValue) const;
    Normal3f FindOneNormal3f(const std::string &name, const Normal3f &defaultValue) const;
    Spectrum FindOneSpectrum(const std::string &name, const Spectrum &defaultValue) const;
    std::string FindOneString(const std::string &name, const std::string &defaultValue) const;
    std::string FindOneFilename(const std::string &name, const std::string &defaultValue) const;
    // No default value for textures.
    std::string FindTexture(const std::string &name) const;

    // Return all the values stored by the item of the given name.
    const bool *FindBool(const std::string &name, int *n) const;
    const int *FindInt(const std::string &name, int *n) const;
    const Float *FindFloat(const std::string &name, int *n) const;
    const Point2f *FindPoint2f(const std::string &name, int *n) const;
    const Vector2f *FindVector2f(const std::string &name, int *n) const;
    const Point3f *FindPoint3f(const std::string &name, int *n) const;
    const Vector3f *FindVector3f(const std::string &name, int *n) const;
    const Normal3f *FindNormal3f(const std::string &name, int *n) const;
    const Spectrum *FindSpectrum(const std::string &name, int *n) const;
    const std::string *FindString(const std::string &name, int *n) const;

    // Add an item with all the given input values under the given name.
    void AddBool(const std::string &name, std::unique_ptr<bool[]> values, int nValues);
    void AddInt(const std::string &name, std::unique_ptr<int[]> values, int nValues);
    void AddFloat(const std::string &name, std::unique_ptr<Float[]> values, int nValues);
    void AddPoint2f(const std::string &name, std::unique_ptr<Point2f[]> values, int nValues);
    void AddVector2f(const std::string &name, std::unique_ptr<Vector2f[]> values, int nValues);
    void AddPoint3f(const std::string &name, std::unique_ptr<Point3f[]> values, int nValues);
    void AddVector3f(const std::string &name, std::unique_ptr<Vector3f[]> values, int nValues);
    void AddNormal3f(const std::string &name, std::unique_ptr<Normal3f[]> values, int nValues);
    // AddRGBSpectrum and AddXYZSpectrum take 3 floating-point values.
    void AddRGBSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues);
    void AddXYZSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues);
    // Takes pairs of temperature (Kelvin) and scale. 
    void AddBlackbodySpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues);
    // Takes pairs of wavelength and SPD values at each wavelength.
    void AddSampledSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues);
    // Reads pairs of wavelength and SPD values at each wavelength from files.
    void AddSampledSpectrumFiles(const std::string &name, const char **names, int nValues);
    void AddString(const std::string &name, std::unique_ptr<std::string[]> values, int nValues);
    void AddTexture(const std::string &name, const std::string &value);

    // Remove the item of the given name.
    bool EraseBool(const std::string &name);
    bool EraseInt(const std::string &name);
    bool EraseFloat(const std::string &name);
    bool ErasePoint2f(const std::string &name);
    bool EraseVector2f(const std::string &name);
    bool ErasePoint3f(const std::string &name);
    bool EraseVector3f(const std::string &name);
    bool EraseNormal3f(const std::string &name);
    bool EraseSpectrum(const std::string &name);
    bool EraseString(const std::string &name);
    bool EraseTexture(const std::string &name);

    // Clear all the parameter vectors.
    void Clear();

    // Warn when an item of the parameter set was never looked up.
    void ReportUnused() const;

private:
    std::vector<std::shared_ptr<ParamSetItem<bool>>> bools;
    std::vector<std::shared_ptr<ParamSetItem<int>>> ints;
    std::vector<std::shared_ptr<ParamSetItem<Float>>> floats;
    std::vector<std::shared_ptr<ParamSetItem<Point2f>>> point2fs;
    std::vector<std::shared_ptr<ParamSetItem<Vector2f>>> vector2fs;
    std::vector<std::shared_ptr<ParamSetItem<Point3f>>> point3fs;
    std::vector<std::shared_ptr<ParamSetItem<Vector3f>>> vector3fs;
    std::vector<std::shared_ptr<ParamSetItem<Normal3f>>> normals;
    std::vector<std::shared_ptr<ParamSetItem<Spectrum>>> spectra;
    std::vector<std::shared_ptr<ParamSetItem<std::string>>> strings;
    std::vector<std::shared_ptr<ParamSetItem<std::string>>> textures;
    static std::map<std::string, Spectrum> cachedSpectra;
};

class TextureParams {
private:
    std::map<std::string, std::shared_ptr<Texture<Float>>> &floatTextures;
    std::map<std::string, std::shared_ptr<Texture<Spectrum>>> &spectrumTextures;
    const ParamSet &geomParams;
    const ParamSet &materialParams;

public:
    TextureParams(
        const ParamSet &geomParams, const ParamSet &materialParams,
        std::map<std::string, std::shared_ptr<Texture<Float>>> &fTex,
        std::map<std::string, std::shared_ptr<Texture<Spectrum>>> &sTex
    ) : floatTextures(fTex),
        spectrumTextures(sTex),
        geomParams(geomParams),
        materialParams(materialParams)
    {}

    // The following return the texture of the input name. The search is done in the 
    // following order:
    // 1. The shape (a whole texture).
    // 2. The material (a whole texture).
    // 3. A constant value from the shape (a single value).
    // 4. A constant value from the material (a single value).
    // 5. The input default value (a single value).

    std::shared_ptr<Texture<Spectrum>> GetSpectrumTexture(
        const std::string &n, const Spectrum &default
    ) const;

    std::shared_ptr<Texture<Spectrum>> GetSpectrumTextureOrNull(
        const std::string &n
    ) const;

    std::shared_ptr<Texture<Float>> GetFloatTexture(
        const std::string &n, Float def
    ) const;

    std::shared_ptr<Texture<Float>> GetFloatTextureOrNull(
        const std::string &n
    ) const;

    // The following return the 1st value of the corresponding type associated with 
    // the input name, or the provided default value if not found. The shape is searched
    // first, then the material.
    //
    // These are not texture functions. These are functions that make TextureParams act
    // as a regular ParamSet.

    Float FindFloat(const std::string &n, Float d) const {
        return geomParams.FindOneFloat(n, materialParams.FindOneFloat(n, d));
    }

    std::string FindString(const std::string &n, const std::string &d = "") const {
        return geomParams.FindOneString(n, materialParams.FindOneString(n, d));
    }

    std::string FindFilename(const std::string &n, const std::string &d = "") const {
        return geomParams.FindOneFilename(n, materialParams.FindOneFilename(n, d));
    }

    int FindInt(const std::string &n, int d) const {
        return geomParams.FindOneInt(n, materialParams.FindOneInt(n, d));
    }

    bool FindBool(const std::string &n, bool d) const {
        return geomParams.FindOneBool(n, materialParams.FindOneBool(n, d));
    }

    Point3f FindPoint3f(const std::string &n, const Point3f &d) const {
        return geomParams.FindOnePoint3f(n, materialParams.FindOnePoint3f(n, d));
    }

    Vector3f FindVector3f(const std::string &n, const Vector3f &d) const {
        return geomParams.FindOneVector3f(n, materialParams.FindOneVector3f(n, d));
    }

    Normal3f FindNormal3f(const std::string &n, const Normal3f &d) const {
        return geomParams.FindOneNormal3f(n, materialParams.FindOneNormal3f(n, d));
    }

    Spectrum FindSpectrum(const std::string &n, const Spectrum &d) const {
        return geomParams.FindOneSpectrum(n, materialParams.FindOneSpectrum(n, d));
    }
};

#endif // CPBRT_CORE_PARAMSET_H