#include "core/light.h"
#include "core/spectrum.h"
#include "core/transform.h"

class PointLight : public Light {
private:
    // Light space origin in world coordinates.
    const Point3f pLight;

    // Constant intensity (flux (power) per unit solid angle) because PointLight is isotropic.    
    const Spectrum I;

public:
    // A PointLight's energy is characterized by a positional delta distribution: a given point may
    // receive radiant intensity from this light only from a single direction.
    PointLight(
        const Transform &LightToWorld,
        // TODO: implement MediumInterface.
        const MediumInterface &mediumInterface,
        const Spectrum &I
    ) : Light((int) LightFlags::DeltaPosition, LightToWorld, mediumInterface),
        pLight(LightToWorld(Point3f(0, 0, 0))),
        I(I)
    {}

    // Computes incident radiance at the point of Interaction.
    // 
    // PointLight implements Sample_Li to conform with the other types of lights whose
    // radiometric quantity is indeed radiance and so that light transport algorithms don't
    // have to have special cases. The "radiance" returned by PointLights has an implicit
    // delta distribution: a single direction. So for a given surface point, PointLight doesn't
    // need to be sampled at multiple incident directions.
    //
    // The returned "radiance" is actually intensity divided by the square of the distance between
    // the Interaction and the light. This gives the quantity the same units as radiance, W/(m^2 sr),
    // while also accounting for fall-off.
    Spectrum Sample_Li(
        const Interaction &it,
        const Point2f &u,
        Vector3f *wi,
        Float *pdf,
        VisibilityTester *vis
    ) const;

    // Computes total emitted power (flux) by integrating intensity over the entire sphere of 
    // directions.
    Spectrum Power() const;
};