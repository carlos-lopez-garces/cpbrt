#include "cpbrt.h"
#include "paramset.h"

// Initialization options stored for global access.
Options CpbrtOptions;

enum class APIState {
    // Before cpbrtInit() and after cpbrtCleanup(). No other API calls are legal.
    Uninitialized,
    // Outside a world block: outside cpbrtWorldBegin()/cpbrtWorldEnd().
    OptionsBlock,
    // Inside cpbrtWorldBegin()/cpbrtWorldEnd().
    WorldBlock
};

static APIState currentApiState = APIState::Uninitialized;

void cpbrtInit(const Options &opt);

void cpbrtCleanup();

// Transformations API.

void cpbrtActiveTransformAll();

void cpbrtActiveTransformEndTime();

void cpbrtActiveTransformStartTime();

// Initializes all the current transformation matrices (CTMs) as identity matrices.
void cpbrtIdentity();

void cpbrtTranslate(Float dx, Float dy, Float dz);

// Rotate active current transformation matrices about the given axis (ax, ay, az).
void cpbrtRotate(Float angle, Float ax, Float ay, Float az);

void cpbrtScale(Float sx, Float sy, Float sz);

void cpbrtLookAt(
    // Camera-space origin in world space.
    Float ex, Float ey, Float ez,
    // Look-at point.
    Float lx, Float ly, Float lz,
    // Up vector. 
    Float ux, Float uy, Float uz
);

// Multiply active current transformations by the input matrix.
void cpbrtConcatTransform(Float transform[16]);

// Set active current transformations to the input transformation.
void cpbrtTransform(Float transform[16]);

// Saves the current transformations under the given name.
void cpbrtCoordinateSystem(const std::string &name);

// Sets the named saved transformations as the current ones.
void cpbrtCoordSysTransform(const std::string &name);

// Options API.

void cpbrtPixelFilter(const std::string &name, const ParamSet &params);

void cpbrtFilm(const std::string &type, const ParamSet &params);

void cpbrtSampler(const std::string &name, const ParamSet &params);

void cpbrtAccelerator(const std::string &name, const ParamSet &params);

void cpbrtIntegrator(const std::string &name, const ParamSet &params);

void cpbrtCamera(const std::string &name, const ParamSet &params);