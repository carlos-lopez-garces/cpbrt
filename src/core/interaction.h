#include "cpbrt.h"
#include "geometry.h"

class Interaction {
public:
    Point3f p;

    Normal3f n;

    Float time;

    // An Interaction carries the floating-point error bound of the associated surface or medium.
    Vector3f pError;

    // Stands for 'omega sub o', the outgoing direction of light at the point p. The 0 vector
    // when 'omega sub o' doesn't apply for the Interaction's use case.
    Vector3f wo;

    // Scattering medium present at the point p.
    MediumInterface mediumInterface;

    Interaction(
        const Point3f &p,
        const Normal3f &n,
        const Vector3f &pError,
        const Vector3f &wo,
        Float time,
        const MediumInterface &mediumInterface
    ) : p(p), n(n), time(time), pError(pError), wo(wo), mediumInterface(mediumInterface) {}

    bool isSurfaceInteraction() const {
        return n != Normal3f();
    }
};