#include "geometry.h"
#include "interaction.h"

class Primitive {
public:
    virtual Bounds3f WorldBound() const = 0;

    virtual bool Intersect(const Ray &r, SurfaceInteraction &si) const = 0;
    // P is for "predicate". No intersection details are returned.
    virtual bool IntersectP(const Ray &r) const = 0;

    virtual const AreaLight *GetAreaLight() const = 0;

    virtual const Material *GetMaterial() const = 0;

    // Includes in the SurfaceInteraction the BSDF and/or BSSRDF that describe the 
    // surface at the ray intersection. The MemoryArena allocates memory for the BSDF
    // and BSSRDF.
    virtual void ComputeScatteringFunctions(
        SurfaceInteraction *si,
        MemoryArena &arena,
        TransportMode mode,
        bool allowMultipleLobes
    ) const = 0;
};