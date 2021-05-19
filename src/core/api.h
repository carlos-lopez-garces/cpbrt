#include "cpbrt.h"

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