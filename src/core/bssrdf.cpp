#include "bssrdf.h"
#include "interpolation.h"

// Computes the 1st moment of the Fresnel reflectance function Fr. The ith Fresnel moment
// is an integral involving the Fresnel reflectance function Fr over the hemisphere of 
// directions. It is approximated here using a polynomial of degree 5 (a Taylor polynomial?).
Float FresnelMoment1(Float eta) {
    Float eta2 = eta * eta;
    Float eta3 = eta2 * eta;
    Float eta4 = eta3 * eta;
    Float eta5 = eta4 * eta;

    if (eta < 1) {
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 + 2.49277f * eta4 - 0.68441f * eta5;
    }
    return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 - 1.27198f * eta4 + 0.12746f * eta5;
}

// Computes the 2nd moment of the Fresnel reflectance function Fr.
Float FresnelMoment2(Float eta) {
    Float eta2 = eta * eta;
    Float eta3 = eta2 * eta;
    Float eta4 = eta3 * eta;
    Float eta5 = eta4 * eta;

    if (eta < 1) {
        return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 + 0.07883f * eta4 + 0.04860f * eta5;
    }
    
    Float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
    return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 + 458.843f * r_eta + 404.557f * eta - 189.519f * eta2 + 54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
}

// Computes the radial profile of BSSRDF using spline-based interpolation of the tabulated
// samples of r and rho.
Spectrum TabulatedBSSRDF::Sr(Float r) const {
    Spectrum Sr(0.f);

    // Do spline-based interpolation for each of the channels of the spectrum.
    for (int ch = 0; ch < Spectrum::nSamples; ++ch) {
        // Reduce the dimensionality of Sr(eta, g, rho, sigma_t, r) by fixing sigma_t=1 and
        // bundling sigma_t and r into a unitless optical radius r_optical = sigma_t * r. 
        // This is effectively a change of variable that requires a scaling factor (a Jacobian
        // determinant): Sr(eta, g, rho, sigma_t, r) = sigma_t^2 * Sr(eta, g, rho, 1, r_optical).
        Float rOptical = r * sigma_t[ch];

        // Compute spline weights to interpolate BSSRDF on channel ch. The weights are the
        // coefficients of a polynomial interpolator of a function f (the Sr table in this case):
        //
        // p(x) = w0*f(x_-1) + w1*f(x_0) + w2*f(x_1) + w3*f(x_2).
        //
        // An interpolator is obtained for rho and another for optical radius.

        // An offset is the index (into the row) that corresponds to the start endpoint of the
        // subinterval (a pair of samples or r or rho) that contains the lookup value (input r
        // and TabulatedBSSRDF::rho).  
        int rhoOffset;
        int radiusOffset;
        // Interpolator weights for rho and optical radius.
        Float rhoWeights[4];
        Float radiusWeights[4];
        if (
            !CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch], &rhoOffset, rhoWeights)
            || !CatmullRomWeights(table.mOpticalRadiusSamples, table.opticalRadiusSamples.get(), rOptical, &radiusOffset, radiusWeights)
        ) {
            // rho or r are not in the domain of the tabulated radial profile Sr.
            continue;
        }

        // Set BSSRDF value Sr[ch] using tensor spline interpolation.
        Float sr = 0;
        // 4 weights.
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Float weight = rhoWeights[i] * radiusWeights[j];
                if (weight != 0) {
                    sr += weight * table.EvalProfile(rhoOffset + i, radiusOffset + j);
                }
            }
        }

        // Cancel marginal PDF factor of 2*Pi*r_optical from tabulatd BSSRDF profile. It
        // was introduced when evaluating the tabulated profile and is there for facilitating
        // importance sampling, but it's not part of the definition of Sr.
        if (rOptical != 0) {
            sr /= 2 * Pi * rOptical;
        }

        Sr[ch] = sr;
    }

    // Transform unitless value into world space units. The Sr value thus far computed is unitless
    // as a result of a change of variable r_optical = sigma_t * r that reduces the dimensionality
    // of Sr (see above).
    Sr *= sigma_t * sigma_t;

    return Sr.Clamp();
}