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
};