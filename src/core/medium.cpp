#include "medium.h"
#include "geometry.h"

struct MeasuredSS {
    const char *name;
    Float sigma_prime_s[3];
    Float sigma_a[3];
};

static MeasuredSS SubsurfaceParameterTable[] = {
    {
        "Cream", {7.38, 5.47, 3.15}, {0.0002, 0.0028, 0.0163},
    },
    {
        "Ketchup", {0.18, 0.07, 0.03}, {0.061, 0.97, 1.45},
    },
    {
        "Marble", {2.19, 2.62, 3.00}, {0.0021, 0.0041, 0.0071},
    },
    {
        "Regular Milk", {4.5513, 5.8294, 7.136}, {0.0015333, 0.0046, 0.019933}
    },
    {
        "Espresso", {0.72378, 0.84557, 1.0247}, {4.7984, 6.5751, 8.8493}
    }
};

bool GetMediumScatteringProperties(
    const std::string &name, Spectrum *sigma_a, Spectrum *sigma_prime_s
) {
    for (MeasuredSS &mss : SubsurfaceParameterTable) {
        if (name == mss.name) {
            *sigma_a = Spectrum::FromRGB(mss.sigma_a);
            *sigma_prime_s = Spectrum::FromRGB(mss.sigma_prime_s);
            return true;
        }
    }
    return false;
}

Float HenyeyGreensteinPhaseFunction::p(const Vector3f &wo, const Vector3f &wi) const {
    return PhaseHG(Dot(wo, wi), g);
}

Float HenyeyGreensteinPhaseFunction::Sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const {
    // Compute  for Henyey–Greenstein sample. 
    Float cosTheta;
    if (std::abs(g) < 1e-3) {
        cosTheta = 1 - 2 * u[0];
    } else {
        Float sqrTerm = (1 - g * g) / (1 - g + 2 * g * u[0]);
        cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }

    // Compute direction wi for Henyey-Greenstein sample. 
    Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
    Float phi = 2 * Pi * u[1];
    Vector3f v1, v2;
    CoordinateSystem(wo, &v1, &v2);
    *wi = SphericalDirection(sinTheta, cosTheta, phi, v1, v2, wo);

    return PhaseHG(cosTheta, g);
}