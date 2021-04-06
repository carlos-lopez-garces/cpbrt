#include "fileutil.h"
#include "paramset.h"

bool ParamSet::FindOneBool(const std::string &name, bool defaultValue) const {
    for (const auto &item : bools) {
        if (item->name == name && item->nValues == 1) {
            item->lookedUp = true;
            return item->values[0];
        }
    }
    return defaultValue;
}

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