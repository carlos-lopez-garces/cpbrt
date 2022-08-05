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
    },
    {
        "White Zinfandel", {1.7501e-05, 1.9069e-05, 1.288e-05}, {0.012072, 0.016184, 0.019843}
    },
    {
        "Merlot", {2.1129e-05, 0, 0}, {0.11632, 0.25191, 0.29434}
    },
    {
        "Budweiser Beer", {2.4356e-05, 2.4079e-05, 1.0564e-05}, {0.011492, 0.024911, 0.057786}
    },
    {
        "Coors Light Beer", {5.0922e-05, 4.301e-05, 0}, {0.006164, 0.013984, 0.034983}
    },
    {
        "Clorox", {0.0024035, 0.0031373, 0.003991}, {0.0033542, 0.014892, 0.026297}
    },
    {
        "Skin1", {0.74, 0.88, 1.01}, {0.032, 0.17, 0.48},
    },
    {
        "Skin2", {1.09, 1.59, 1.79}, {0.013, 0.070, 0.145},
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
    // Compute  for Henyey-Greenstein sample. 
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

Float HenyeyGreensteinPhaseFunction::Pdf(Vector3f wo, Vector3f wi) const {
    return p(wo, wi);
}