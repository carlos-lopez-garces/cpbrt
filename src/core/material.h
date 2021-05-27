#ifndef CPBRT_CORE_MATERIAL_H
#define CPBRT_CORE_MATERIAL_H

#include "interaction.h"
#include "memory.h"

// TODO: describe.
enum class TransportMode {
  Radiance,
  Importance
};

class Material {
public:
  virtual void ComputeScatteringFunctions(
      // Differential geometry of surface-ray intersection point.
      SurfaceInteraction *si,
      // For allocating BSDFs.
      MemoryArena &arena,
      TransportMode mode,
      bool allowMultipleLobes
  ) const = 0;  
};

#endif // CPBRT_CORE_MATERIAL_H