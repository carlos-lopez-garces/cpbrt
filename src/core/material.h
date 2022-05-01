#ifndef CPBRT_CORE_MATERIAL_H
#define CPBRT_CORE_MATERIAL_H

#include "cpbrt.h"
#include "memory.h"

// TODO: describe.
// Did the ray start at the camera or at a light source?
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