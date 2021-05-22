#include "api.h"
#include "error.h"
#include "spectrum.h"
#include "transform.h"

// API local classes.

constexpr int MaxTransforms = 2;
constexpr int StartTransformBits = 1 << 0;
constexpr int EndTransformBits   = 1 << 1;
constexpr int AllTransformsBits = (1 << MaxTransforms) - 1;

// Stores an array of transformations.
struct TransformSet {
private:
    Transform t[MaxTransforms];

public:
    Transform &operator[](int i) {
        return t[i];
    }

    // Returns the inverses of the transformations in a TansformSet. 
    friend TansformSet Inverse(const TransformSet &ts) {
        TransformSet tInv;
        for (int i = 0; i < MaxTransforms; ++i) {
            tInv.t[i] = Inverse(ts.t[i]);
        }
        return tInv;
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
};

// Stack of attributes.
struct GraphicsState {

};

// API static data.

// Current transformation matrices (CTMs).
static TransformSet curTransform;
static uint32_t activeTransformBits = AllTransformsBits;

// Saved CTMs.
static std::map<std::string, TransformSet> namedCoordinateSystems;

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
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookat;);
}

void cpbrtConcatTransform(Float transform[16]) {
    VERIFY_INITIALIZED("ConcatTransform");
    FOR_ACTIVE_TRANSFORMS(
        curTransform[i] =
            curTransform[i] *
            Transform(Matrix4x4(tr[0], tr[4], tr[8], tr[12], tr[1], tr[5],
                                tr[9], tr[13], tr[2], tr[6], tr[10], tr[14],
                                tr[3], tr[7], tr[11], tr[15])););
}

void cpbrtTransform(Float transform[16]) {
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