#include "bssrdf.h"

// Computes the 1st moment of the Fresnel reflectance function Fr. The ith Fresnel moment
// is an integral involving the Fresnel reflectance function Fr over the hemisphere of 
// directions. It is approximated here using a polynomial of degree 5 (a Taylor polynomial?).
Float FresnelMoment1(Float reciprocalEta) {
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
Float FresnelMoment2(Float reciprocalEta) {
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