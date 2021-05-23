#include <map>

#include "fileutil.h"
#include "paramset.h"
#include "spectrum.h"
#include "textures/constant.h"

// ParamSet data members.

std::map<std::string, Spectrum> ParamSet::cachedSpectra;

// ParamSet macros.

#define ADD_PARAM_TYPE(T, vec) \
    (vec).emplace_back(new ParamSetItem<T>(name, std::move(values), nValues));

bool ParamSet::FindOneBool(const std::string &name, bool defaultValue) const {
    for (const auto &item : bools) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

// ParamSet methods.

int ParamSet::FindOneInt(const std::string &name, int defaultValue) const {
    for (const auto &item : ints) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Float ParamSet::FindOneFloat(const std::string &name, Float defaultValue) const {
    for (const auto &item : floats) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Point2f ParamSet::FindOnePoint2f(const std::string &name, const Point2f &defaultValue) const {
    for (const auto &item : point2fs) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Vector2f ParamSet::FindOneVector2f(const std::string &name, const Vector2f &defaultValue) const {
    for (const auto &item : vector2fs) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Point3f ParamSet::FindOnePoint3f(const std::string &name, const Point3f &defaultValue) const {
    for (const auto &item : point3fs) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Vector3f ParamSet::FindOneVector3f(const std::string &name, const Vector3f &defaultValue) const {
    for (const auto &item : vector3fs) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Normal3f ParamSet::FindOneNormal3f(const std::string &name, const Normal3f &defaultValue) const {
    for (const auto &item : normals) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

Spectrum ParamSet::FindOneSpectrum(const std::string &name, const Spectrum &defaultValue) const {
    for (const auto &item : spectra) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

std::string ParamSet::FindOneString(const std::string &name, const std::string &defaultValue) const {
    for (const auto &item : strings) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

std::string ParamSet::FindOneFilename(const std::string &name, const std::string &defaultValue) const {
    std::string filename = FindOneString(name, "");
    if (filename == "") {
        return defaultValue;
    }
    filename = AbsolutePath(ResolveFilename(filename));
    return filename;
}

std::string ParamSet::FindTexture(const std::string &name) const {
    std::string defaultValue = "";
    for (const auto &item : textures) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

const bool *ParamSet::FindBool(const std::string &name, int *n) const {
    for (const auto &item : bools) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const int *ParamSet::FindInt(const std::string &name, int *n) const {
    for (const auto &item : ints) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Float *ParamSet::FindFloat(const std::string &name, int *n) const {
    for (const auto &item : floats) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Point2f *ParamSet::FindPoint2f(const std::string &name, int *n) const {
    for (const auto &item : point2fs) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Vector2f *ParamSet::FindVector2f(const std::string &name, int *n) const {
    for (const auto &item : vector2fs) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Point3f *ParamSet::FindPoint3f(const std::string &name, int *n) const {
    for (const auto &item : point3fs) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Vector3f *ParamSet::FindVector3f(const std::string &name, int *n) const {
    for (const auto &item : vector3fs) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Normal3f *ParamSet::FindNormal3f(const std::string &name, int *n) const {
    for (const auto &item : normals) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const Spectrum *ParamSet::FindSpectrum(const std::string &name, int *n) const {
    for (const auto &item : spectra) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

const std::string *ParamSet::FindString(const std::string &name, int *n) const {
    for (const auto &item : strings) {
        if (item->name == name) {
            *n = item->nValues;
            item->lookedUp = true;
            return item->values.get();
        }
    }
    return nullptr;
}

void ParamSet::AddBool(const std::string &name, std::unique_ptr<bool[]> values, int nValues) {
    EraseBool(name);
    ADD_PARAM_TYPE(bool, bools);
}

void ParamSet::AddInt(const std::string &name, std::unique_ptr<int[]> values, int nValues) {
    EraseInt(name);
    ADD_PARAM_TYPE(int, ints);
}

void ParamSet::AddFloat(const std::string &name, std::unique_ptr<Float[]> values, int nValues) {
    EraseFloat(name);
    floats.emplace_back(new ParamSetItem<Float>(name, std::move(values), nValues));
}

void ParamSet::AddPoint2f(const std::string &name, std::unique_ptr<Point2f[]> values, int nValues) {
    ErasePoint2f(name);
    ADD_PARAM_TYPE(Point2f, point2fs);
}

void ParamSet::AddVector2f(const std::string &name, std::unique_ptr<Vector2f[]> values, int nValues) {
    EraseVector2f(name);
    ADD_PARAM_TYPE(Vector2f, vector2fs);
}

void ParamSet::AddPoint3f(const std::string &name, std::unique_ptr<Point3f[]> values, int nValues) {
    ErasePoint3f(name);
    ADD_PARAM_TYPE(Point3f, point3fs);
}

void ParamSet::AddVector3f(const std::string &name, std::unique_ptr<Vector3f[]> values, int nValues) {
    EraseVector3f(name);
    ADD_PARAM_TYPE(Vector3f, vector3fs);
}

void ParamSet::AddNormal3f(const std::string &name, std::unique_ptr<Normal3f[]> values, int nValues) {
    EraseNormal3f(name);
    ADD_PARAM_TYPE(Normal3f, normals);
}

void ParamSet::AddRGBSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues) {
    EraseSpectrum(name);
    nValues /= 3;
    std::unique_ptr<Spectrum[]> s(new Spectrum[nValues]);
    for (int i = 0; i < nValues; ++i) s[i] = Spectrum::FromRGB(&values[3 * i]);
    std::shared_ptr<ParamSetItem<Spectrum>> psi(new ParamSetItem<Spectrum>(name, std::move(s), nValues));
    spectra.push_back(psi);
}

void ParamSet::AddXYZSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues) {
    EraseSpectrum(name);
    nValues /= 3;
    std::unique_ptr<Spectrum[]> s(new Spectrum[nValues]);
    for (int i = 0; i < nValues; ++i) s[i] = Spectrum::FromXYZ(&values[3 * i]);
    std::shared_ptr<ParamSetItem<Spectrum>> psi(new ParamSetItem<Spectrum>(name, std::move(s), nValues));
    spectra.push_back(psi);
}

void ParamSet::AddBlackbodySpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues) {
    EraseSpectrum(name);
    nValues /= 2;
    std::unique_ptr<Spectrum[]> s(new Spectrum[nValues]);
    std::unique_ptr<Float[]> v(new Float[nCIESamples]);
    for (int i = 0; i < nValues; ++i) {
        BlackbodyNormalized(CIE_lambda, nCIESamples, values[2 * i], v.get());
        s[i] = values[2 * i + 1] * Spectrum::FromSampled(CIE_lambda, v.get(), nCIESamples);
    }
    std::shared_ptr<ParamSetItem<Spectrum>> psi(
        new ParamSetItem<Spectrum>(name, std::move(s), nValues));
    spectra.push_back(psi);
}

void ParamSet::AddSampledSpectrum(const std::string &name, std::unique_ptr<Float[]> values, int nValues) {
    EraseSpectrum(name);
    nValues /= 2;
    std::unique_ptr<Float[]> wl(new Float[nValues]);
    std::unique_ptr<Float[]> v(new Float[nValues]);
    for (int i = 0; i < nValues; ++i) {
        wl[i] = values[2 * i];
        v[i] = values[2 * i + 1];
    }
    std::unique_ptr<Spectrum[]> s(new Spectrum[1]);
    s[0] = Spectrum::FromSampled(wl.get(), v.get(), nValues);
    std::shared_ptr<ParamSetItem<Spectrum>> psi(new ParamSetItem<Spectrum>(name, std::move(s), 1));
    spectra.push_back(psi);
}

void ParamSet::AddSampledSpectrumFiles(const std::string &name, const char **names, int nValues) {
    EraseSpectrum(name);
    std::unique_ptr<Spectrum[]> s(new Spectrum[nValues]);
    for (int i = 0; i < nValues; ++i) {
        std::string fn = AbsolutePath(ResolveFilename(names[i]));
        if (cachedSpectra.find(fn) != cachedSpectra.end()) {
            s[i] = cachedSpectra[fn];
            continue;
        }

        std::vector<Float> vals;
        if (!ReadFloatFile(fn.c_str(), &vals)) {
            Warning(
                "Unable to read SPD file \"%s\".  Using black distribution.",
                fn.c_str());
            s[i] = Spectrum(0.);
        } else {
            if (vals.size() % 2) {
                Warning(
                    "Extra value found in spectrum file \"%s\". "
                    "Ignoring it.",
                    fn.c_str());
            }
            std::vector<Float> wls, v;
            for (size_t j = 0; j < vals.size() / 2; ++j) {
                wls.push_back(vals[2 * j]);
                v.push_back(vals[2 * j + 1]);
            }
            s[i] = Spectrum::FromSampled(&wls[0], &v[0], wls.size());
        }
        cachedSpectra[fn] = s[i];
    }

    std::shared_ptr<ParamSetItem<Spectrum>> psi(new ParamSetItem<Spectrum>(name, std::move(s), nValues));
    spectra.push_back(psi);
}

void ParamSet::AddString(const std::string &name, std::unique_ptr<std::string[]> values, int nValues) {
    EraseString(name);
    ADD_PARAM_TYPE(std::string, strings);
}

void ParamSet::AddTexture(const std::string &name, const std::string &value) {
    EraseTexture(name);
    std::unique_ptr<std::string[]> str(new std::string[1]);
    str[0] = value;
    std::shared_ptr<ParamSetItem<std::string>> psi(new ParamSetItem<std::string>(name, std::move(str), 1));
    textures.push_back(psi);
}

bool ParamSet::EraseBool(const std::string &name) {
    for (size_t i = 0; i < bools.size(); ++i) {
        if (bools[i]->name == name) {
            bools.erase(bools.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseInt(const std::string &name) {
    for (size_t i = 0; i < ints.size(); ++i) {
        if (ints[i]->name == name) {
            ints.erase(ints.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseFloat(const std::string &name) {
    for (size_t i = 0; i < floats.size(); ++i) {
        if (floats[i]->name == name) {
            floats.erase(floats.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::ErasePoint2f(const std::string &name) {
    for (size_t i = 0; i < point2fs.size(); ++i) {
        if (point2fs[i]->name == name) {
            point2fs.erase(point2fs.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseVector2f(const std::string &name) {
    for (size_t i = 0; i < vector2fs.size(); ++i) {
        if (vector2fs[i]->name == name) {
            vector2fs.erase(vector2fs.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::ErasePoint3f(const std::string &name) {
    for (size_t i = 0; i < point3fs.size(); ++i) {
        if (point3fs[i]->name == name) {
            point3fs.erase(point3fs.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseVector3f(const std::string &name) {
    for (size_t i = 0; i < vector3fs.size(); ++i) {
        if (vector3fs[i]->name == name) {
            vector3fs.erase(vector3fs.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseNormal3f(const std::string &name) {
    for (size_t i = 0; i < normals.size(); ++i) {
        if (normals[i]->name == name) {
            normals.erase(normals.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseSpectrum(const std::string &name) {
    for (size_t i = 0; i < spectra.size(); ++i) {
        if (spectra[i]->name == name) {
            spectra.erase(spectra.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseString(const std::string &name) {
    for (size_t i = 0; i < strings.size(); ++i) {
        if (strings[i]->name == name) {
            strings.erase(strings.begin() + i);
            return true;
        }
    }
    return false;
}

bool ParamSet::EraseTexture(const std::string &name) {
    for (size_t i = 0; i < textures.size(); ++i) {
        if (textures[i]->name == name) {
            textures.erase(textures.begin() + i);
            return true;
        }
    }
    return false;
}

void ParamSet::Clear() {
#define DEL_PARAMS(name) (name).erase((name).begin(), (name).end())
    DEL_PARAMS(ints);
    DEL_PARAMS(bools);
    DEL_PARAMS(floats);
    DEL_PARAMS(point2fs);
    DEL_PARAMS(vector2fs);
    DEL_PARAMS(point3fs);
    DEL_PARAMS(vector3fs);
    DEL_PARAMS(normals);
    DEL_PARAMS(spectra);
    DEL_PARAMS(strings);
    DEL_PARAMS(textures);
#undef DEL_PARAMS
}

void ParamSet::ReportUnused() const {
#define CHECK_UNUSED(v)                                                 \
    for (size_t i = 0; i < (v).size(); ++i)                             \
        if (!(v)[i]->lookedUp)                                          \
            Warning("Parameter \"%s\" not used", (v)[i]->name.c_str())
    CHECK_UNUSED(ints);
    CHECK_UNUSED(bools);
    CHECK_UNUSED(floats);
    CHECK_UNUSED(point2fs);
    CHECK_UNUSED(vector2fs);
    CHECK_UNUSED(point3fs);
    CHECK_UNUSED(vector3fs);
    CHECK_UNUSED(normals);
    CHECK_UNUSED(spectra);
    CHECK_UNUSED(strings);
    CHECK_UNUSED(textures);
}

std::shared_ptr<Texture<Spectrum>> TextureParams::GetSpectrumTexture(
    const std::string &n,
    const Spectrum &def
) const {
    std::shared_ptr<Texture<Spectrum>> tex = GetSpectrumTextureOrNull(n);
    if (tex)
        return tex;
    else
        return std::make_shared<ConstantTexture<Spectrum>>(def);
}

std::shared_ptr<Texture<Spectrum>> TextureParams::GetSpectrumTextureOrNull(
    const std::string &n
) const {
    // Check the shape parameters first.
    std::string name = geomParams.FindTexture(n);
    if (name.empty()) {
        int count;
        const Spectrum *s = geomParams.FindSpectrum(n, &count);
        if (s) {
            if (count > 1)
                Warning("Ignoring excess values provided with parameter \"%s\"",
                        n.c_str());
            return std::make_shared<ConstantTexture<Spectrum>>(*s);
        }

        name = materialParams.FindTexture(n);
        if (name.empty()) {
            int count;
            const Spectrum *s = materialParams.FindSpectrum(n, &count);
            if (s) {
                if (count > 1)
                    Warning("Ignoring excess values provided with parameter \"%s\"",
                            n.c_str());
                return std::make_shared<ConstantTexture<Spectrum>>(*s);
            }
        }

        if (name.empty())
            return nullptr;
    }

    // We have a texture name, from either the shape or the material's
    // parameters.
    if (spectrumTextures.find(name) != spectrumTextures.end())
        return spectrumTextures[name];
    else {
        Error("Couldn't find spectrum texture named \"%s\" for parameter \"%s\"",
              name.c_str(), n.c_str());
        return nullptr;
    }
}

std::shared_ptr<Texture<Float>> TextureParams::GetFloatTexture(
    const std::string &n,
    Float def
) const {
    std::shared_ptr<Texture<Float>> tex = GetFloatTextureOrNull(n);
    if (tex)
        return tex;
    else
        return std::make_shared<ConstantTexture<Float>>(def);
}

std::shared_ptr<Texture<Float>> TextureParams::GetFloatTextureOrNull(
    const std::string &n
) const {
    // Check the shape parameters first.
    std::string name = geomParams.FindTexture(n);
    if (name.empty()) {
        int count;
        const Float *s = geomParams.FindFloat(n, &count);
        if (s) {
            if (count > 1)
                Warning("Ignoring excess values provided with parameter \"%s\"",
                        n.c_str());
            return std::make_shared<ConstantTexture<Float>>(*s);
        }

        name = materialParams.FindTexture(n);
        if (name.empty()) {
            int count;
            const Float *s = materialParams.FindFloat(n, &count);
            if (s) {
                if (count > 1)
                    Warning("Ignoring excess values provided with parameter \"%s\"",
                            n.c_str());
                return std::make_shared<ConstantTexture<Float>>(*s);
            }
        }

        if (name.empty())
            return nullptr;
    }

    // We have a texture name, from either the shape or the material's
    // parameters.
    if (floatTextures.find(name) != floatTextures.end())
        return floatTextures[name];
    else {
        Error("Couldn't find float texture named \"%s\" for parameter \"%s\"",
              name.c_str(), n.c_str());
        return nullptr;
    }
}