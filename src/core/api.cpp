#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "film.h"
#include "medium.h"

#include "accelerators/bvh.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/directlighting.h"
#include "integrators/path.h"
#include "lights/point.h"
#include "materials/matte.h"
#include "materials/plastic.h"
#include "samplers/stratified.h"
#include "shapes/sphere.h"
#include "textures/constant.h"

// API local classes.

constexpr int MaxTransforms = 2;
constexpr int StartTransformBits = 1 << 0;
constexpr int EndTransformBits   = 1 << 1;
constexpr int AllTransformsBits = (1 << MaxTransforms) - 1;

// Static functions.

std::shared_ptr<Texture<Float>> MakeFloatTexture(
    const std::string &name,
    const Transform &tex2world,
    const TextureParams &tp
) {
    Texture<Float> *tex = nullptr;
    
    // TODO: add other types.
    if (name == "constant")
        tex = CreateConstantFloatTexture(tex2world, tp);
    else
        Warning("Float texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();

    return std::shared_ptr<Texture<Float>>(tex);
}

std::shared_ptr<Texture<Spectrum>> MakeSpectrumTexture(
    const std::string &name,
    const Transform &tex2world,
    const TextureParams &tp
) {
    Texture<Spectrum> *tex = nullptr;

    // TODO: add other types.
    if (name == "constant")
        tex = CreateConstantSpectrumTexture(tex2world, tp);
    else
        Warning("Spectrum texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();

    return std::shared_ptr<Texture<Spectrum>>(tex);
}

std::shared_ptr<Material> MakeMaterial(const std::string &name, const TextureParams &mp) {
    Material *material = nullptr;

    // TODO: add other types.
    if (name == "" || name == "none")
        return nullptr;
    else if (name == "matte")
        material = CreateMatteMaterial(mp);
    else if (name == "plastic")
        material = CreatePlasticMaterial(mp);
    else {
        Warning("Material \"%s\" unknown. Using \"matte\".", name.c_str());
        material = CreateMatteMaterial(mp);
    }
    mp.ReportUnused();

    if (!material) Error("Unable to create material \"%s\"", name.c_str());
    return std::shared_ptr<Material>(material);
}

std::shared_ptr<Light> MakeLight(
    const std::string &name,
    const ParamSet &paramSet,
    const Transform &light2world,
    const MediumInterface &mediumInterface
) {
    std::shared_ptr<Light> light;

    // TODO: add other types.
    if (name == "point") {
        light =
            CreatePointLight(light2world, mediumInterface.outside, paramSet);
    }
    else {
        Warning("Light \"%s\" unknown.", name.c_str()); 
    }
    paramSet.ReportUnused();
    
    return light;
}

std::vector<std::shared_ptr<Shape>> MakeShapes(
    const std::string &name,
    const Transform *ObjectToWorld,
    const Transform *WorldToObject,
    bool reverseOrientation,
    const ParamSet &paramSet
) {
    std::vector<std::shared_ptr<Shape>> shapes;
    std::shared_ptr<Shape> s;

    // TODO: add other types.
    if (name == "sphere") {
        s = CreateSphereShape(ObjectToWorld, WorldToObject, reverseOrientation, paramSet);
    }
    
    if (s != nullptr) {
        shapes.push_back(s);
    } else {
        Warning("Shape \"%s\" unknown.", name.c_str());
    }

    return shapes;
}

// Stores an array of transformations.
struct TransformSet {
private:
    Transform t[MaxTransforms];

public:
    Transform &operator[](int i) {
        return t[i];
    }

    // Returns the inverses of the transformations in a TansformSet. 
    friend TransformSet Inverse(const TransformSet &ts) {
        TransformSet tInv;
        for (int i = 0; i < MaxTransforms; ++i) {
            tInv.t[i] = Inverse(ts.t[i]);
        }
        return tInv;
    }

    bool IsAnimated() const {
        for (int i = 0; i < MaxTransforms - 1; ++i) {
            if (t[i] != t[i+1]) {
                return true;
            }
        }
        return false;
    }
};

// Caches the inverse of a transformation.
class TransformCache {
private:
    std::map<Transform, std::pair<Transform *, Transform *>> cache;
    MemoryArena arena;

public:
    void Lookup(const Transform &t, Transform **tCached, Transform **tCachedInverse) {
        auto iter = cache.find(t);
        if (iter == cache.end()) {
            // Cache miss. Allocate and store.
            Transform *tr = arena.Alloc<Transform>();
            *tr = t;
            Transform *tInv = arena.Alloc<Transform>();
            *tInv = Transform(Inverse(t));
            cache[t] = std::make_pair(tr, tInv);
            iter  = cache.find(t);
        }

        if (tCached) {
            *tCached = iter->second.first;
        }
        if (tCachedInverse) {
            *tCachedInverse = iter->second.second;
        }
    }

    void Clear() {
        arena.Reset();
        cache.erase(cache.begin(), cache.end());
    }
};

// Options that can be set in the OptionsBlock state.
struct RenderOptions {
    std::string FilterName = "box";
    ParamSet FilterParams;

    std::string FilmName = "image";
    ParamSet FilmParams;

    std::string SamplerName = "stratified";
    ParamSet SamplerParams;

    std::string AcceleratorName = "bvh";
    ParamSet AcceleratorParams;

    std::string IntegratorName = "path";
    ParamSet IntegratorParams;

    std::string CameraName = "perspective";
    ParamSet CameraParams;
    TransformSet CameraToWorld;

    std::vector<std::shared_ptr<Light>> lights;

    // TODO: describe.
    std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
};

// An instance of a material.
struct MaterialInstance {
    MaterialInstance() = default;

    MaterialInstance(
        const std::string &name,
        const std::shared_ptr<Material> &mtl,
        ParamSet params
    ) : name(name), material(mtl), params(std::move(params))
    {}

    std::string name;
    std::shared_ptr<Material> material;
    ParamSet params;
};

// Stack of attributes.
struct GraphicsState {
    using FloatTextureMap = std::map<std::string, std::shared_ptr<Texture<Float>>>;
    std::shared_ptr<FloatTextureMap> floatTextures;
    bool floatTexturesShared = false;

    using SpectrumTextureMap = std::map<std::string, std::shared_ptr<Texture<Spectrum>>>;
    std::shared_ptr<SpectrumTextureMap> spectrumTextures;
    bool spectrumTexturesShared = false;

    std::shared_ptr<MaterialInstance> currentMaterial;
    ParamSet materialParams;
    std::string material = "matte";

    using NamedMaterialMap = std::map<std::string, std::shared_ptr<MaterialInstance>>;
    std::shared_ptr<NamedMaterialMap> namedMaterials;
    bool namedMaterialsShared = false;

    ParamSet areaLightParams;
    std::string areaLight;

    // TODO: describe.
    bool reverseOrientation = false; 

    GraphicsState() 
        : floatTextures(std::make_shared<FloatTextureMap>()),
          spectrumTextures(std::make_shared<SpectrumTextureMap>()),
          namedMaterials(std::make_shared<NamedMaterialMap>()) {
        
        ParamSet empty;
        TextureParams tp(empty, empty, *floatTextures, *spectrumTextures);
        std::shared_ptr<Material> mtl(CreateMatteMaterial(tp));
        currentMaterial = std::make_shared<MaterialInstance>("matte", mtl, ParamSet());
    }

    MediumInterface CreateMediumInterface();

    std::shared_ptr<Material> GetMaterialForShape(const ParamSet &params);
};

MediumInterface GraphicsState::CreateMediumInterface() {
    // TODO: implement.
    MediumInterface m;
    return m;
}

// TODO: explain.
bool shapeMaySetMaterialParameters(const ParamSet &ps) {
    for (const auto &param : ps.textures)
        if (param->name != "alpha" && param->name != "shadowalpha")
            return true;
    for (const auto &param : ps.floats)
        if (param->nValues == 1 && param->name != "radius")
            return true;
    for (const auto &param : ps.strings)
        if (param->nValues == 1 && param->name != "filename" &&
            param->name != "type" && param->name != "scheme")
            return true;
    for (const auto &param : ps.bools)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.ints)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.normals)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.spectra)
        if (param->nValues == 1)
            return true;

    return false;
}

std::shared_ptr<Material> GraphicsState::GetMaterialForShape(const ParamSet &params) {
    if (shapeMaySetMaterialParameters(params)) {
        TextureParams mp(params, currentMaterial->params, *floatTextures, *spectrumTextures);
        return MakeMaterial(currentMaterial->name, mp);
    } else {
        return currentMaterial->material;
    }
}

// API static data.

// Current transformation matrices (CTMs).
static TransformSet curTransform;
static uint32_t activeTransformBits = AllTransformsBits;

// Saved CTMs.
static std::map<std::string, TransformSet> namedCoordinateSystems;

// Cached inverses of transformations.
static TransformCache transformCache;

// Scene-wide global options set in the OptionsBlock state.
static std::unique_ptr<RenderOptions> renderOptions;

static GraphicsState graphicsState;
static std::vector<GraphicsState> pushedGraphicsStates;
static std::vector<TransformSet> pushedTransforms;
static std::vector<uint32_t> pushedActiveTransformBits;

// API macros.

// API functions that are only legal while in the initialized state call this macro.
#define VERIFY_INITIALIZED(func)                                   \
    if (!(CpbrtOptions.cat || CpbrtOptions.toPly) &&               \
        currentApiState == APIState::Uninitialized) {              \
        Error(                                                     \
            "cpbrtInit() must be called before calling \"%s()\". " \
            "Ignoring.",                                           \
            func);                                                 \
        return;                                                    \
    } else /* swallow trailing semicolon */

// API functions that are only legal in an options block call this macro.
#define VERIFY_OPTIONS(func)                             \
    VERIFY_INITIALIZED(func);                            \
    if (!(CpbrtOptions.cat || CpbrtOptions.toPly) &&     \
        currentApiState == APIState::WorldBlock) {       \
        Error(                                           \
            "Options cannot be set inside world block; " \
            "\"%s\" not allowed.  Ignoring.",            \
            func);                                       \
        return;                                          \
    } else /* swallow trailing semicolon */

// API functions that are only legal in a world block call this macro.
#define VERIFY_WORLD(func)                                   \
    VERIFY_INITIALIZED(func);                                \
    if (!(CpbrtOptions.cat || CpbrtOptions.toPly) &&         \
        currentApiState == APIState::OptionsBlock) {         \
        Error(                                               \
            "Scene description must be inside world block; " \
            "\"%s\" not allowed. Ignoring.",                 \
            func);                                           \
        return;                                              \
    } else /* swallow trailing semicolon */

// Executes the expression on active transformations only.
#define FOR_ACTIVE_TRANSFORMS(expr)                   \
    for (int i = 0; i < MaxTransforms; ++i)           \
        if (activeTransformBits & (1 << i)) { expr }

void cpbrtInit(const Options &opt) {
    CpbrtOptions = opt;

    if (currentApiState != APIState::Uninitialized) {
        Error("cpbrtInit() has already been called.");
    }
    currentApiState = APIState::OptionsBlock;
    renderOptions.reset(new RenderOptions());

    graphicsState = GraphicsState();

    SampledSpectrum::Init();
}

void cpbrtCleanup() {
    if (currentApiState == APIState::Uninitialized) {
        Error("cpbrtCleanup() called without cpbrtInit().");
    } else if (currentApiState == APIState::WorldBlock) {
        Error("cpbrtCleanup() called while inside world block.");
    }
    currentApiState = APIState::Uninitialized;
    renderOptions.reset(nullptr);
}

// Transformations API.

void cpbrtActiveTransformAll() {
    activeTransformBits = AllTransformsBits;
}

void cpbrtActiveTransformEndTime() {
    activeTransformBits = EndTransformBits;
}

void cpbrtActiveTransformStartTime() {
    activeTransformBits = StartTransformBits;
}

void cpbrtIdentity() {
    VERIFY_INITIALIZED("Identity");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform(););
}

void cpbrtTranslate(Float dx, Float dy, Float dz) {
    VERIFY_INITIALIZED("Translate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Translate(Vector3f(dx, dy, dz)););
}

void cpbrtRotate(Float angle, Float ax, Float ay, Float az) {
    VERIFY_INITIALIZED("Rotate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Rotate(angle, Vector3f(ax, ay, az)););
}

void cpbrtScale(Float sx, Float sy, Float sz) {
    VERIFY_INITIALIZED("Scale");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Scale(sx, sy, sz););
}

void cpbrtLookAt(
    Float ex, Float ey, Float ez, 
    Float lx, Float ly, Float lz, 
    Float ux, Float uy, Float uz
) {
    VERIFY_INITIALIZED("LookAt");
    Transform lookAt = LookAt(Point3f(ex, ey, ez), Point3f(lx, ly, lz), Vector3f(ux, uy, uz));
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookAt;);
}

void cpbrtConcatTransform(Float tr[16]) {
    VERIFY_INITIALIZED("ConcatTransform");
    FOR_ACTIVE_TRANSFORMS(
        curTransform[i] =
            curTransform[i] *
            Transform(Matrix4x4(tr[0], tr[4], tr[8], tr[12], tr[1], tr[5],
                                tr[9], tr[13], tr[2], tr[6], tr[10], tr[14],
                                tr[3], tr[7], tr[11], tr[15])););
}

void cpbrtTransform(Float tr[16]) {
    VERIFY_INITIALIZED("Transform");
    FOR_ACTIVE_TRANSFORMS(
        curTransform[i] = Transform(Matrix4x4(
            tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2],
            tr[6], tr[10], tr[14], tr[3], tr[7], tr[11], tr[15])););
}

void cpbrtCoordinateSystem(const std::string &name) {
    VERIFY_INITIALIZED("CoordinateSystem");
    namedCoordinateSystems[name] = curTransform;
}

void cpbrtCoordSysTransform(const std::string &name) {
    VERIFY_INITIALIZED("CoordSysTransform");
    if (namedCoordinateSystems.find(name) != namedCoordinateSystems.end()) {
        curTransform = namedCoordinateSystems[name];
    } else {
        Warning("Couldn't find named coordinate system \"%s\"", name.c_str());
    }
}

// Options API.

void cpbrtPixelFilter(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("PixelFilter");
    renderOptions->FilterName = name;
    renderOptions->FilterParams = params;
}

void cpbrtFilm(const std::string &type, const ParamSet &params) {
    VERIFY_OPTIONS("Film");
    renderOptions->FilmName = type;
    renderOptions->FilmParams = params;
}

void cpbrtSampler(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Sampler");
    renderOptions->SamplerName = name;
    renderOptions->SamplerParams = params;
}

void cpbrtAccelerator(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Accelerator");
    renderOptions->AcceleratorName = name;
    renderOptions->AcceleratorParams = params;
}

void cpbrtIntegrator(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Integrator");
    renderOptions->IntegratorName = name;
    renderOptions->IntegratorParams = params;
}

void cpbrtCamera(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Camera");
    renderOptions->CameraName = name;
    renderOptions->CameraParams = params;
    renderOptions->CameraToWorld = Inverse(curTransform);
    namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
}

// World API.

void cpbrtWorldBegin() {
    VERIFY_OPTIONS("WorldBegin");
    currentApiState = APIState::WorldBlock;
    for (int i = 0; i < MaxTransforms; ++i) {
        // Reset all CTMs as identity matrices.
        curTransform[i] = Transform();
    }
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems["world"] = curTransform;
}

void pbrtObjectBegin(const std::string &name) {
    VERIFY_WORLD("ObjectBegin");
    cpbrtAttributeBegin();
    if (renderOptions->currentInstance)
        Error("ObjectBegin called inside of instance definition");
    renderOptions->instances[name] = std::vector<std::shared_ptr<Primitive>>();
    renderOptions->currentInstance = &renderOptions->instances[name];
}

void pbrtObjectEnd() {
    VERIFY_WORLD("ObjectEnd");
    if (!renderOptions->currentInstance)
        Error("ObjectEnd called outside of instance definition");
    renderOptions->currentInstance = nullptr;
    cpbrtAttributeEnd();
}

void cpbrtAttributeBegin() {
    VERIFY_WORLD("AttributeBegin");
    pushedGraphicsStates.push_back(graphicsState);
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}

void cpbrtAttributeEnd() {
    VERIFY_WORLD("AttributeEnd");
    if (!pushedGraphicsStates.size()) {
        Error("Unmatched cpbrtAttributeEnd() encountered. Ignoring it."); 
        return;
    }
    graphicsState = pushedGraphicsStates.back();
    pushedGraphicsStates.pop_back();
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}

void cpbrtTransformBegin() {
    VERIFY_WORLD("TransformBegin");
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}

void cpbrtTransformEnd() {
    VERIFY_WORLD("TransformEnd");
    if (!pushedTransforms.size()) {
        Error(
            "Unmatched cpbrtTransformEnd() encountered. Ignoring it.");
        return;
    }
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}

void cpbrtReverseOrientation() {
    VERIFY_WORLD("ReverseOrientation");
    graphicsState.reverseOrientation = !graphicsState.reverseOrientation;
}

void cpbrtTexture(
    const std::string &name,
    const std::string &type,
    const std::string &texname,
    const ParamSet &params
) {
    VERIFY_WORLD("Texture");

    TextureParams tp(params, params, *graphicsState.floatTextures, *graphicsState.spectrumTextures);

    if (type == "float") {
        if (graphicsState.floatTextures->find(name) != graphicsState.floatTextures->end()) {
            Warning("Texture \"%s\" being redefined", name.c_str());
        }
            
        std::shared_ptr<Texture<Float>> ft = MakeFloatTexture(texname, curTransform[0], tp);
        if (ft) {
            if (graphicsState.floatTexturesShared) {
                graphicsState.floatTextures =
                    std::make_shared<GraphicsState::FloatTextureMap>(*graphicsState.floatTextures);
                graphicsState.floatTexturesShared = false;
            }
            (*graphicsState.floatTextures)[name] = ft;
        }
    } else if (type == "color" || type == "spectrum") {
        if (graphicsState.spectrumTextures->find(name) != graphicsState.spectrumTextures->end()) {
            Warning("Texture \"%s\" being redefined", name.c_str());
        }

        std::shared_ptr<Texture<Spectrum>> st = MakeSpectrumTexture(texname, curTransform[0], tp);
        if (st) {
            if (graphicsState.spectrumTexturesShared) {
                graphicsState.spectrumTextures =
                    std::make_shared<GraphicsState::SpectrumTextureMap>(*graphicsState.spectrumTextures);
                graphicsState.spectrumTexturesShared = false;
            }
            (*graphicsState.spectrumTextures)[name] = st;
        }
    } else
        Error("Texture type \"%s\" unknown.", type.c_str());
}

void cpbrtMaterial(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Material");
    ParamSet emptyParams;
    TextureParams mp(params, emptyParams, *graphicsState.floatTextures, *graphicsState.spectrumTextures);
    std::shared_ptr<Material> mtl = MakeMaterial(name, mp);
    graphicsState.currentMaterial = std::make_shared<MaterialInstance>(name, mtl, params);
}

void cpbrtMakeNamedMaterial(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("MakeNamedMaterial");
    
    ParamSet emptyParams;
    TextureParams mp(params, emptyParams, *graphicsState.floatTextures, *graphicsState.spectrumTextures);
    std::string matName = mp.FindString("type");
    
    if (matName == "") {
        Error("No parameter string \"type\" found in MakeNamedMaterial");
    }

    std::shared_ptr<Material> mtl = MakeMaterial(matName, mp);
    
    if (graphicsState.namedMaterials->find(name) != graphicsState.namedMaterials->end()) {
        Warning("Named material \"%s\" redefined.", name.c_str());
    }

    if (graphicsState.namedMaterialsShared) {
        graphicsState.namedMaterials =
            std::make_shared<GraphicsState::NamedMaterialMap>(*graphicsState.namedMaterials);
        graphicsState.namedMaterialsShared = false;
    }
    
    (*graphicsState.namedMaterials)[name] = std::make_shared<MaterialInstance>(matName, mtl, params);
}

void cpbrtNamedMaterial(const std::string &name) {
    VERIFY_WORLD("NamedMaterial");

    auto iter = graphicsState.namedMaterials->find(name);
    if (iter == graphicsState.namedMaterials->end()) {
        Error("NamedMaterial \"%s\" unknown.", name.c_str());
        return;
    }
    graphicsState.currentMaterial = iter->second;
}

void cpbrtLightSource(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("LightSource");

     MediumInterface mi = graphicsState.CreateMediumInterface();

    std::shared_ptr<Light> lt = MakeLight(name, params, curTransform[0], mi);

    if (!lt) {
        Error("LightSource: light type \"%s\" unknown.", name.c_str());
    } else {
        renderOptions->lights.push_back(lt);
    }
}

void cpbrtAreaLightSource(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("AreaLightSource");
    graphicsState.areaLight = name;
    graphicsState.areaLightParams = params;
    // Unlike cpbrtLightSource, the area light is not created yet: the shapes
    // that make up its geometry need to be created first.
}

void cpbrtShape(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Shape");
    std::vector<std::shared_ptr<Primitive>> prims;
    std::vector<std::shared_ptr<AreaLight>> areaLights;
    if (!curTransform.IsAnimated()) {
        // Initialize prims and areaLights for static shape.

        // Create shapes.
        Transform *ObjToWorld;
        Transform *WorldToObj;
        transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
        std::vector<std::shared_ptr<Shape>> shapes 
            = MakeShapes(name, ObjToWorld, WorldToObj, graphicsState.reverseOrientation, params);
        if (shapes.size() == 0) {
            return;
        }

        std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
        params.ReportUnused();
        
        MediumInterface mi = graphicsState.CreateMediumInterface();

        prims.reserve(shapes.size());

        for (auto s : shapes) {
            // Possibly create area light.
            std::shared_ptr<AreaLight> areaLight;
            if (graphicsState.areaLight != "") {
                // TODO: implement.
            }

            prims.push_back(std::make_shared<GeometricPrimitive>(s, mtl, areaLight, mi));
        }
    } else {
        // TODO: implement. Also TransformedPrimitive and AnimatedTransforms.
    }
    // TODO: Add prims and areaLights to scene or current instance.
}