#include "api.h"
#include "error.h"
#include "spectrum.h"

// API macros.

// API functions that are only legal while in the initialized state call this macro.
#define VERIFY_INITIALIZED(func)                                  \
    if (!(CpbrtOptions.cat || CpbrtOptions.toPly) &&              \
        currentApiState == APIState::Uninitialized) {             \
        Error(                                                    \
            "pbrtInit() must be called before calling \"%s()\". " \
            "Ignoring.",                                          \
            func);                                                \
        return;                                                   \
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

void cpbrtInit(const Options &opt) {
    CpbrtOptions = opt;

    if (currentApiState != APIState::Uninitialized) {
        Error("cpbrtInit() has already been called.");
    }
    currentApiState = APIState::OptionsBlock;

    SampledSpectrum::Init();
}

void cpbrtCleanup() {
    if (currentApiState == APIState::Uninitialized) {
        Error("cpbrtCleanup() called without cpbrtInit().");
    } else if (currentApiState == APIState::WorldBlock) {
        Error("cpbrtCleanup() called while inside world block.");
    }
    currentApiState = APIState::Uninitialized;
}