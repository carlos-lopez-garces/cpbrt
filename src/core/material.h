#include "interaction.h"
#include "memory.h"

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