#include "cpbrt.h"
#include "microfacet.h"
#include "reflection.h"

// TODO.
Float MicrofacetDistribution::Pdf(const Vector3f &wo,
                                  const Vector3f &wh) const {
    if (sampleVisibleArea)
        return D(wh) * G1(wo) * AbsDot(wo, wh) / AbsCosTheta(wo);
    else
        return D(wh) * AbsCosTheta(wh);
}

// TODO.
static void BeckmannSample11(Float cosThetaI, Float U1, Float U2,
                             Float *slope_x, Float *slope_y) {
    /* Special case (normal incidence) */
    if (cosThetaI > .9999) {
        Float r = std::sqrt(-std::log(1.0f - U1));
        Float sinPhi = std::sin(2 * Pi * U2);
        Float cosPhi = std::cos(2 * Pi * U2);
        *slope_x = r * cosPhi;
        *slope_y = r * sinPhi;
        return;
    }

    /* The original inversion routine from the paper contained
       discontinuities, which causes issues for QMC integration
       and techniques like Kelemen-style MLT. The following code
       performs a numerical inversion with better behavior */
    Float sinThetaI =
        std::sqrt(std::max((Float)0, (Float)1 - cosThetaI * cosThetaI));
    Float tanThetaI = sinThetaI / cosThetaI;
    Float cotThetaI = 1 / tanThetaI;

    /* Search interval -- everything is parameterized
       in the Erf() domain */
    Float a = -1, c = Erf(cotThetaI);
    Float sample_x = std::max(U1, (Float)1e-6f);

    /* Start with a good initial guess */
    // Float b = (1-sample_x) * a + sample_x * c;

    /* We can do better (inverse of an approximation computed in
     * Mathematica) */
    Float thetaI = std::acos(cosThetaI);
    Float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
    Float b = c - (1 + c) * std::pow(1 - sample_x, fit);

    /* Normalization factor for the CDF */
    static const Float SQRT_PI_INV = 1.f / std::sqrt(Pi);
    Float normalization =
        1 /
        (1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));

    int it = 0;
    while (++it < 10) {
        /* Bisection criterion -- the oddly-looking
           Boolean expression are intentional to check
           for NaNs at little additional cost */
        if (!(b >= a && b <= c)) b = 0.5f * (a + c);

        /* Evaluate the CDF and its derivative
           (i.e. the density function) */
        Float invErf = ErfInv(b);
        Float value =
            normalization *
                (1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) -
            sample_x;
        Float derivative = normalization * (1 - invErf * tanThetaI);

        if (std::abs(value) < 1e-5f) break;

        /* Update bisection intervals */
        if (value > 0)
            c = b;
        else
            a = b;

        b -= value / derivative;
    }

    /* Now convert back into a slope value */
    *slope_x = ErfInv(b);

    /* Simulate Y component */
    *slope_y = ErfInv(2.0f * std::max(U2, (Float)1e-6f) - 1.0f);
}

// TODO.
static Vector3f BeckmannSample(const Vector3f &wi, Float alpha_x, Float alpha_y,
                               Float U1, Float U2) {
    // 1. stretch wi
    Vector3f wiStretched =
        Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    Float slope_x, slope_y;
    BeckmannSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

    // 3. rotate
    Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    return Normalize(Vector3f(-slope_x, -slope_y, 1.f));
}

Float BeckmannDistribution::D(const Vector3f &wh) const {
    Float tan2Theta = Tan2Theta(wh);

    if (std::isinf(tan2Theta)) {
        // wh is a perfectly grazing angle. The Beckmann-Spizzichino
        // distribution function approaches 0 as tan2Theta approaches infinity.
        return 0.f;
    }

    Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);

    // Obtain the RMS slope (alpha) that corresponds to normal wh by doing an ellipse
    // interpolation of the RMS slopes (alphaX and alphaY) that correspond to the 
    // x-axis-perpendicular and y-axis-perpendicular azimuthal orientation angles phiX and phiY.
    // The interpolator is the actual azimuthal orientation angle of the wh normal.
    Float ellipseInterpolatedAlpha = (Cos2Phi(wh) / (alphaX*alphaX)) + (Sin2Phi(wh) / (alphaY*alphaY)); 

    return std::exp(-tan2Theta * ellipseInterpolatedAlpha) / (Pi * alphaX * alphaY * cos4Theta);
}

Float BeckmannDistribution::Lambda(const Vector3f &w) const {
    Float absTanTheta = std::abs(TanTheta(w));

    if (std::isinf(absTanTheta)) {
        // Perfectly grazing angle.
        return 0.f;
    };

    Float interpolatedAlpha = std::sqrt(Cos2Phi(w)*alphaX*alphaX + Sin2Phi(w)*alphaY*alphaY);

    // Polynomial approximation to Lambda.
    Float a = 1 / (interpolatedAlpha * absTanTheta);
    if (a >= 1.6f) {
        return 0.f;
    }
    return (1 - 1.259f*a + 0.396f*a*a) / (3.535f*a + 2.181f*a*a);
}

// TODO.
Vector3f BeckmannDistribution::Sample_wh(const Vector3f &wo,
                                         const Point2f &u) const {
    if (!sampleVisibleArea) {
        // Sample full distribution of normals for Beckmann distribution

        // Compute $\tan^2 \theta$ and $\phi$ for Beckmann distribution sample
        Float tan2Theta, phi;
        if (alphaX == alphaY) {
            Float logSample = std::log(1 - u[0]);
            tan2Theta = -alphaX * alphaX * logSample;
            phi = u[1] * 2 * Pi;
        } else {
            // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
            // distribution
            Float logSample = std::log(1 - u[0]);

            phi = std::atan(alphaY / alphaX *
                            std::tan(2 * Pi * u[1] + 0.5f * Pi));
            if (u[1] > 0.5f) phi += Pi;
            Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            Float alphaX2 = alphaX * alphaX, alphaY2 = alphaY * alphaY;
            tan2Theta = -logSample /
                        (cosPhi * cosPhi / alphaX2 + sinPhi * sinPhi / alphaY2);
        }

        // Map sampled Beckmann angles to normal direction _wh_
        Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
        Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
        Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
        return wh;
    } else {
        // Sample visible area of normals for Beckmann distribution
        Vector3f wh;
        bool flip = wo.z < 0;
        wh = BeckmannSample(flip ? -wo : wo, alphaX, alphaY, u[0], u[1]);
        if (flip) wh = -wh;
        return wh;
    }
}