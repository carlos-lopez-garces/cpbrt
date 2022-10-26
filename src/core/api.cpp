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
#include "cameras/realistic.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/directlighting.h"
#include "integrators/path.h"
#include "integrators/volpath.h"
#include "lights/diffuse.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "materials/coateddiffuse.h"
#include "materials/glass.h"
#include "materials/matte.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/plastic.h"
#include "materials/subsurface.h"
#include "materials/ward.h"
#include "media/grid.h"
#include "media/homogeneous.h"
#include "samplers/stratified.h"
#include "shapes/sphere.h"
#include "shapes/triangle.h"
#include "textures/checkerboard.h"
#include "textures/constant.h"

// Initialization options stored for global access.
Options CpbrtOptions;

constexpr int MaxTransforms = 2;
constexpr int StartTransformBits = 1 << 0;
constexpr int EndTransformBits   = 1 << 1;
constexpr int AllTransformsBits = (1 << MaxTransforms) - 1;

// API local classes.

// Stores an array of transformations.
struct TransformSet {
private:
    Transform t[MaxTransforms];

public:
    Transform &operator[](int i) {
        return t[i];
    }

    const Transform &operator[](int i) const {
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
  public:
    TransformCache()
        : hashTable(512), hashTableOccupancy(0) {}

    Transform *Lookup(const Transform &t) {
        int offset = Hash(t) & (hashTable.size() - 1);
        int step = 1;
        while (true) {
            if (!hashTable[offset] || *hashTable[offset] == t)
                break;
            offset = (offset + step * step) & (hashTable.size() - 1);
            ++step;
        }
        Transform *tCached = hashTable[offset];
        if (!tCached) {
            tCached = arena.Alloc<Transform>();
            *tCached = t;
            Insert(tCached);
        }
        return tCached;
    }

    void Clear() {
        hashTable.clear();
        hashTable.resize(512);
        hashTableOccupancy = 0;
        arena.Reset();
    }

  private:
    void Insert(Transform *tNew);
    void Grow();

    static uint64_t Hash(const Transform &t) {
        const char *ptr = (const char *)(&t.GetMatrix());
        size_t size = sizeof(Matrix4x4);
        uint64_t hash = 14695981039346656037ull;
        while (size > 0) {
            hash ^= *ptr;
            hash *= 1099511628211ull;
            ++ptr;
            --size;
        }
        return hash;
    }

    std::vector<Transform *> hashTable;
    int hashTableOccupancy;
    MemoryArena arena;
};

void TransformCache::Insert(Transform *tNew) {
    if (++hashTableOccupancy == hashTable.size() / 2)
        Grow();

    int baseOffset = Hash(*tNew) & (hashTable.size() - 1);
    for (int nProbes = 0;; ++nProbes) {
        // Quadratic probing.
        int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
        if (hashTable[offset] == nullptr) {
            hashTable[offset] = tNew;
            return;
        }
    }
}

void TransformCache::Grow() {
    std::vector<Transform *> newTable(2 * hashTable.size());

    for (Transform *tEntry : hashTable) {
        if (!tEntry) continue;

        int baseOffset = Hash(*tEntry) & (hashTable.size() - 1);
        for (int nProbes = 0;; ++nProbes) {
            int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
            if (newTable[offset] == nullptr) {
                newTable[offset] = tEntry;
                break;
            }
        }
    }

    std::swap(hashTable, newTable);
}

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

    std::string currentInsideMedium;
    std::string currentOutsideMedium;

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

    std::vector<std::shared_ptr<Primitive>> primitives;

    std::map<std::string, std::shared_ptr<Medium>> namedMedia;
    bool haveScatteringMedia = false;

    // TODO: describe.
    std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;

    Integrator *MakeIntegrator() const;

    Scene *MakeScene();

    Camera *MakeCamera() const;
};

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

// Static functions.

Camera *MakeCamera(
    const std::string &name,
    const ParamSet &paramSet,
    const TransformSet &cam2worldSet,
    Film *film
) {
    Camera *camera = nullptr;
    MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
    static_assert(MaxTransforms == 2, "TransformCache assumes only two transforms");
    Transform *cam2world[2] = {
        transformCache.Lookup(cam2worldSet[0]),
        transformCache.Lookup(cam2worldSet[1])
    };

    Transform cam2World(cam2world[0]->GetMatrix());
    if (name == "perspective") {
        camera = CreatePerspectiveCamera(paramSet, cam2World, film, mediumInterface.outside);
    } else if (name == "orthographic") {
        camera = CreateOrthographicCamera(paramSet, cam2World, film, mediumInterface.outside);
    } else if (name == "realistic") {
        camera = CreateRealisticCamera(paramSet, cam2World, film, mediumInterface.outside);
    } else {
        Warning("Camera \"%s\" unknown.", name.c_str());
    }
    paramSet.ReportUnused();
    
    return camera;
}

std::shared_ptr<Sampler> MakeSampler(
    const std::string &name,
    const ParamSet &paramSet,
    const Film *film
) {
    Sampler *sampler = nullptr;

    // TODO: add other types.
    if (name == "stratified") {
        sampler = CreateStratifiedSampler(paramSet);
    }
    else {
        Warning("Sampler \"%s\" unknown.", name.c_str());
    }
    
    paramSet.ReportUnused();
    return std::shared_ptr<Sampler>(sampler);
}

std::unique_ptr<Filter> MakeFilter(const std::string &name,
                                   const ParamSet &paramSet) {
    Filter *filter = nullptr;

    if (name == "box")
        filter = CreateBoxFilter(paramSet);
    else if (name == "gaussian")
        filter = CreateGaussianFilter(paramSet);
    else if (name == "mitchell")
        filter = CreateMitchellFilter(paramSet);
    else if (name == "sinc")
        filter = CreateSincFilter(paramSet);
    else if (name == "triangle")
        filter = CreateTriangleFilter(paramSet);
    else {
        Error("Filter \"%s\" unknown.", name.c_str());
        exit(1);
    }
    paramSet.ReportUnused();

    return std::unique_ptr<Filter>(filter);
}

Film *MakeFilm(
    const std::string &name,
    const ParamSet &paramSet,
    std::unique_ptr<Filter> filter
) {
    Film *film = nullptr;
    if (name == "image") {
        film = CreateFilm(paramSet, std::move(filter));
    } else {
        Warning("Film \"%s\" unknown.", name.c_str());
    }
    paramSet.ReportUnused();

    return film;
}

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
    if (name == "constant") {
        tex = CreateConstantSpectrumTexture(tex2world, tp);
    } else if (name == "checkerboard") {
        tex = CreateCheckerboardSpectrumTexture(tex2world, tp); 
    }
    else {
        Warning("Spectrum texture \"%s\" unknown.", name.c_str());
    }
    tp.ReportUnused();

    return std::shared_ptr<Texture<Spectrum>>(tex);
}

std::shared_ptr<Material> MakeMaterial(const std::string &name, const TextureParams &mp) {
    Material *material = nullptr;

    // TODO: add other types.
    if (name == "" || name == "none") {
        return nullptr;
    } else if (name == "matte") {
        material = CreateMatteMaterial(mp);
    } else if (name == "mirror") {
        material = CreateMirrorMaterial(mp);
    } else if (name == "plastic") {
        material = CreatePlasticMaterial(mp);
    } else if (name == "metal" || name == "conductor") {
        material = CreateMetalMaterial(mp);
    } else if (name == "glass") {
        material = CreateGlassMaterial(mp); 
    } else if (name == "coateddiffuse") {
        material = CreateCoatedDiffuseMaterial(mp);
    } else if (name == "subsurface") {
        material = CreateSubsurfaceMaterial(mp);
    } else if (name == "ward") {
        material = CreateWardMaterial(mp);
    } else {
        Warning("Material \"%s\" unknown. Using \"matte\".", name.c_str());
        material = CreateMatteMaterial(mp);
    }
    mp.ReportUnused();

    if (!material) Error("Unable to create material \"%s\"", name.c_str());
    return std::shared_ptr<Material>(material);
}

std::shared_ptr<Medium> MakeMedium(
    const std::string &name,
    const ParamSet &paramSet,
    const Transform &medium2world
) {
    Float sig_a_rgb[3] = {.0011f, .0024f, .014f};
    Float sig_s_rgb[3] = {2.55f, 3.21f, 3.77f};

    Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb);
    Spectrum sig_s = Spectrum::FromRGB(sig_s_rgb);

    std::string preset = paramSet.FindOneString("preset", "");
    
    bool found = GetMediumScatteringProperties(preset, &sig_a, &sig_s);
    if (preset != "" && !found) {
        Warning("Material preset \"%s\" not found.  Using defaults.", preset.c_str());
    }

    Float scale = paramSet.FindOneFloat("scale", 1.f);
    Float g = paramSet.FindOneFloat("g", 0.0f);
    sig_a = paramSet.FindOneSpectrum("sigma_a", sig_a) * scale;
    sig_s = paramSet.FindOneSpectrum("sigma_s", sig_s) * scale;

    Medium *m = NULL;
    if (name == "homogeneous") {
        m = new HomogeneousMedium(sig_a, sig_s, g);
    } else if (name == "heterogeneous") {
        int nitems;
        const Float *data = paramSet.FindFloat("density", &nitems);
        if (!data) {
            Error("No \"density\" values provided for heterogeneous medium?");
            return NULL;
        }

        int nx = paramSet.FindOneInt("nx", 1);
        int ny = paramSet.FindOneInt("ny", 1);
        int nz = paramSet.FindOneInt("nz", 1);
        Point3f p0 = paramSet.FindOnePoint3f("p0", Point3f(0.f, 0.f, 0.f));
        Point3f p1 = paramSet.FindOnePoint3f("p1", Point3f(1.f, 1.f, 1.f));

        if (nitems != nx * ny * nz) {
            Error(
                "GridDensityMedium has %d density values; expected nx*ny*nz = "
                "%d",
                nitems, nx * ny * nz);
            return NULL;
        }

        Transform data2Medium = Translate(Vector3f(p0)) * Scale(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
        m = new GridDensityMedium(sig_a, sig_s, g, nx, ny, nz, medium2world * data2Medium, data);
    } else {
        Warning("Medium \"%s\" unknown.", name.c_str());
    }
    paramSet.ReportUnused();
    return std::shared_ptr<Medium>(m);
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
        light = CreatePointLight(light2world, mediumInterface.outside, paramSet);
    } else if (name == "infinite") {
        light = CreateInfiniteLight(light2world, paramSet);
    } else {
        Warning("Light \"%s\" unknown.", name.c_str()); 
    }
    paramSet.ReportUnused();
    
    return light;
}

std::shared_ptr<AreaLight> MakeAreaLight(
    const std::string &name,
    const Transform &LightToWorld,
    const MediumInterface &mediumInterface,
    const ParamSet &params,
    const std::shared_ptr<Shape> &shape
) {
    std::shared_ptr<AreaLight> areaLight;

    // "area" and "diffuse" are synonymous because only diffuse area lights area supported.
    if (name == "area" || name == "diffuse") {
        areaLight = CreateDiffuseAreaLight(LightToWorld, mediumInterface.outside, params, shape);
    } else {
        Warning("Area light \"%s\" unknown.", name.c_str()); 
    }
    params.ReportUnused();

    return areaLight;
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
    } else if (name == "trianglemesh") {
        shapes = CreateTriangleMeshShape(
            ObjectToWorld, WorldToObject, reverseOrientation, paramSet, &*graphicsState.floatTextures
        );
    } else if (name == "plymesh") {
        shapes = CreatePLYMesh(
            ObjectToWorld, WorldToObject, reverseOrientation, paramSet, &*graphicsState.floatTextures
        );
    } else {
        Warning("Shape \"%s\" unknown.", name.c_str());
    }

    return shapes;
}

std::shared_ptr<Primitive> MakeAccelerator(
    const std::string &name,
    std::vector<std::shared_ptr<Primitive>> prims,
    const ParamSet &paramSet) {
    std::shared_ptr<Primitive> accel;

    // TODO: add other types.
    if (name == "bvh") {
        accel = CreateBVHAccelerator(std::move(prims), paramSet);
    }
    else {
        Warning("Accelerator \"%s\" unknown.", name.c_str());
    }

    paramSet.ReportUnused();
    return accel;
}

// RenderOptions function definitions.

Integrator *RenderOptions::MakeIntegrator() const {
    std::shared_ptr<const Camera> camera(MakeCamera());
    if (!camera) {
        Error("Unable to create camera");
        return nullptr;
    }

    std::shared_ptr<Sampler> sampler = MakeSampler(SamplerName, SamplerParams, camera->film);
    if (!sampler) {
        Error("Unable to create sampler.");
        return nullptr;
    }

    Integrator *integrator = nullptr;
    if (IntegratorName == "directlighting") {
        integrator = CreateDirectLightingIntegrator(IntegratorParams, sampler, camera);
    }
    else if (IntegratorName == "path") {
        integrator = CreatePathIntegrator(IntegratorParams, sampler, camera);
    } else if (IntegratorName == "volpath") {
        integrator = CreateVolPathIntegrator(IntegratorParams, sampler, camera);
    } else {
        Error("Integrator \"%s\" unknown.", IntegratorName.c_str());
        return nullptr;
    }

    IntegratorParams.ReportUnused();
    if (lights.empty()) {
        Warning(
            "No light sources defined in scene; "
            "rendering a black image.");
    }

    return integrator;
}

Scene *RenderOptions::MakeScene() {
    std::shared_ptr<Primitive> accelerator
        = ::MakeAccelerator(AcceleratorName, std::move(primitives), AcceleratorParams);
    if (!accelerator) accelerator = std::make_shared<BVHAccel>(primitives);
    Scene *scene = new Scene(accelerator, lights);
    primitives.clear();
    lights.clear();
    return scene;
}

Camera *RenderOptions::MakeCamera() const {
    std::unique_ptr<Filter> filter = MakeFilter(FilterName, FilterParams);
    Film *film = MakeFilm(FilmName, FilmParams, std::move(filter));
    if (!film) {
        Error("Unable to create film.");
        return nullptr;
    }
    Camera *camera = ::MakeCamera(
        CameraName, CameraParams, CameraToWorld, film
    );

    return camera;
}

// GraphicsState function definitions.

MediumInterface GraphicsState::CreateMediumInterface() {
    MediumInterface m;
    if (currentInsideMedium != "") {
        if (renderOptions->namedMedia.find(currentInsideMedium) 
            != renderOptions->namedMedia.end()
        ) {
            m.inside = renderOptions->namedMedia[currentInsideMedium].get();
        }
        else {
            Error("Named medium \"%s\" undefined.", currentInsideMedium.c_str());
        }
    }
    if (currentOutsideMedium != "") {
        if (renderOptions->namedMedia.find(currentOutsideMedium)
            != renderOptions->namedMedia.end()) {
            m.outside = renderOptions->namedMedia[currentOutsideMedium].get();
        } else {
            Error("Named medium \"%s\" undefined.", currentOutsideMedium.c_str());
        }
    }
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

void cpbrtWorldEnd() {
    VERIFY_WORLD("WorldEnd");
    while (pushedGraphicsStates.size()) {
        Warning("Missing end to pbrtAttributeBegin()");
        pushedGraphicsStates.pop_back();
        pushedTransforms.pop_back();
    }
    while (pushedTransforms.size()) {
        Warning("Missing end to pbrtTransformBegin()");
        pushedTransforms.pop_back();
    }

    std::unique_ptr<Integrator> integrator(renderOptions->MakeIntegrator());
    std::unique_ptr<Scene> scene(renderOptions->MakeScene());

    if (scene && integrator) integrator->Render(*scene);

    graphicsState = GraphicsState();
    transformCache.Clear();
    currentApiState = APIState::OptionsBlock;
    renderOptions.reset(new RenderOptions);

    for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems.erase(namedCoordinateSystems.begin(), namedCoordinateSystems.end());
}

void cpbrtObjectBegin(const std::string &name) {
    VERIFY_WORLD("ObjectBegin");
    cpbrtAttributeBegin();
    if (renderOptions->currentInstance)
        Error("ObjectBegin called inside of instance definition");
    renderOptions->instances[name] = std::vector<std::shared_ptr<Primitive>>();
    renderOptions->currentInstance = &renderOptions->instances[name];
}

void cpbrtObjectEnd() {
    VERIFY_WORLD("ObjectEnd");
    if (!renderOptions->currentInstance)
        Error("ObjectEnd called outside of instance definition");
    renderOptions->currentInstance = nullptr;
    cpbrtAttributeEnd();
}

void cpbrtObjectInstance(const std::string &name) {
    VERIFY_WORLD("ObjectInstance");

    if (renderOptions->currentInstance) {
        Error("ObjectInstance can't be called inside instance definition");
        return;
    }
    if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
        Error("Unable to find instance named \"%s\"", name.c_str());
        return;
    }

    std::vector<std::shared_ptr<Primitive>> &in = renderOptions->instances[name];
    if (in.empty()) return;
    if (in.size() > 1) {
        std::shared_ptr<Primitive> accel(
            MakeAccelerator(renderOptions->AcceleratorName, std::move(in), renderOptions->AcceleratorParams)
        );
        if (!accel) accel = std::make_shared<BVHAccel>(in);
        in.clear();
        in.push_back(accel);
    }

    static_assert(MaxTransforms == 2, "TransformCache assumes only two transforms");

    Transform *InstanceToWorld[2] = {
        transformCache.Lookup(curTransform[0]),
        transformCache.Lookup(curTransform[1])
    };

    Transform instanceToWorld(InstanceToWorld[0]->GetMatrix());

    std::shared_ptr<Primitive> prim(std::make_shared<TransformedPrimitive>(in[0], instanceToWorld));

    renderOptions->primitives.push_back(prim);
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

void cpbrtMakeNamedMedium(const std::string &name, const ParamSet &params) {
    VERIFY_INITIALIZED("MakeNamedMedium");

    std::string type = params.FindOneString("type", "");
    if (type == "") {
        Error("No parameter string \"type\" found in MakeNamedMedium");
    } else {
        std::shared_ptr<Medium> medium = MakeMedium(type, params, curTransform[0]);
        if (medium) renderOptions->namedMedia[name] = medium;
    }
}

void cpbrtMediumInterface(const std::string &insideName, const std::string &outsideName) {
    VERIFY_INITIALIZED("MediumInterface");
    graphicsState.currentInsideMedium = insideName;
    graphicsState.currentOutsideMedium = outsideName;
    renderOptions->haveScatteringMedia = true;
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
        Transform *ObjToWorld = transformCache.Lookup(curTransform[0]);
        Transform *WorldToObj = transformCache.Lookup(Inverse(curTransform[0]));
        std::vector<std::shared_ptr<Shape>> shapes 
            = MakeShapes(name, ObjToWorld, WorldToObj, graphicsState.reverseOrientation, params);
        if (shapes.size() == 0) {
            return;
        }

        std::shared_ptr<Material> mtl = graphicsState.GetMaterialForShape(params);
        params.ReportUnused();
        
        MediumInterface mi = graphicsState.CreateMediumInterface();

        prims.reserve(shapes.size());

        for (auto shape : shapes) {
            // Possibly create area light.
            std::shared_ptr<AreaLight> areaLight;
            if (graphicsState.areaLight != "") {
                areaLight = MakeAreaLight(
                    graphicsState.areaLight,
                    curTransform[0],
                    mi,
                    graphicsState.areaLightParams,
                    shape
                );
                if (areaLight) areaLights.push_back(areaLight);
            }

            prims.push_back(std::make_shared<GeometricPrimitive>(shape, mtl, areaLight, mi));
        }
    } else {
        // TODO: implement. Also TransformedPrimitive and AnimatedTransforms.
    }

    if (renderOptions->currentInstance) {
      // TODO: implement object instancing.
    }
    else {
      renderOptions->primitives.insert(renderOptions->primitives.end(), prims.begin(), prims.end());
      
      // Possibly add area light.
      if (areaLights.size()) {
        renderOptions->lights.insert(renderOptions->lights.end(), areaLights.begin(), areaLights.end());
      }
    }
}