#ifndef CPBRT_CORE_API_H
#define CPBRT_CORE_API_H

#include "cpbrt.h"

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

// World API.

void cpbrtWorldBegin();

// TODO: explain.
void cpbrtWorldEnd();

// TODO: explain.
void cpbrtObjectBegin(const std::string &name);

// TODO: explain.
void cpbrtObjectEnd();

// TODO: explain.
void cpbrtObjectInstance(const std::string &name);

// Pushes the current graphics state and transformations onto the stack.
void cpbrtAttributeBegin();

// Pops the current graphics state and transformations from the stack.
void cpbrtAttributeEnd();

// Pushes the CTMs onto the stack.
void cpbrtTransformBegin();

// Pops transformations from the stack.
void cpbrtTransformEnd();

void cpbrtReverseOrientation();

// Records a named texture in the graphics state.
void cpbrtTexture(
    const std::string &name,
    const std::string &type,
    const std::string &texname,
    const ParamSet &params
);

// Records a material in the graphics state.
void cpbrtMaterial(const std::string &name, const ParamSet &params);

// Creates a named material.
void cpbrtMakeNamedMaterial(const std::string &name, const ParamSet &params);

// Records the material of the input name as the current material. cpbrtMakeNamedMaterial
// must have been called to create such named material.
void cpbrtNamedMaterial(const std::string &name);

// Creates a named medium.
void cpbrtMakeNamedMedium(const std::string &name, const ParamSet &params);

void cpbrtMediumInterface(const std::string &insideName, const std::string &outsideName);

void cpbrtLightSource(const std::string &name, const ParamSet &params);

void cpbrtAreaLightSource(const std::string &name, const ParamSet &params);

void cpbrtShape(const std::string &name, const ParamSet &params);

// Main program API.

void cpbrtParseFile(std::string filename);

#endif // CPBRT_CORE_API_H