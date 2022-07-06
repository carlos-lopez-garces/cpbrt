#ifndef CPBRT_CORE_SPECTRUM_H
#define CPBRT_CORE_SPECTRUM_H

#include <cmath>

#include "cpbrt.h"

static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;

enum class SpectrumType { 
    Reflectance,
    Illuminant
};

// Values of the CIE XYZ color-matching curves X(lambda), Y(lambda), Z(lambda)
// for 471 wavelength samples in the [360nm, 830nm] interval.
static const int nCIESamples = 471;
extern const Float CIE_lambda[nCIESamples];
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
// Definite integral of the CIE Y (luminance) color-matching curve over the [360nm, 830nm]
// integration interval.
static const Float CIE_Y_integral = 106.856895;

// Precomputed SPD samples for RGB primaries and colors that are mixes of primaries exclusively,
// used for mapping RGB tristimulus values to SPDs.
//
// Among the infinitely many metameric SPDs that map to the same RGB tristimulus value, the
// sampled SPDs chosen for these colors are smooth (as opposed to spiky) and constant
// when r=g=b.
//
// SampledSpectrums resample these SPDs at their number of samples.
static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
// For reflected light.
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];
// For emitted light. The CIE illumninant D65 SPD is used for white, and the sampled SPDs for
// the rest of the precomputed colors are defined in reference to this white. 
extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];

// Determines whether the samples are sorted in order of increasing wavelength.
extern bool SpectrumSamplesSorted(const Float *lambda, int n);

// Sorts the samples in order of increasing wavelength.
extern void SortSpectrumSamples(Float *lambda, Float *values, int n);

// Computes the average value of the samples contained in the wavelength subinterval
// [lambdaStart, lambdaEnd].
extern Float AverageSpectrumSamples(
    const Float *lambda, 
    const Float *values,
    int n,
    Float lambdaStart,
    Float lambdaEnd
);

extern Float InterpolateSpectrumSamples(const Float *lambda, const Float *values, int n, Float l);

// Computes exitant radiance Le for each of the n wavelengths at temperature T according
// to Planck's law of blackbody emission.
extern void Blackbody(const Float *lambda, int n, Float T, Float *Le);

// Computes normalized exitant radiance for blackbodies.
extern void BlackbodyNormalized(const Float *lambda, int n, Float T, Float *Le);

// Linearly maps XYZ SPD coefficients to RGB SPD coefficients.
inline void XYZToRGB(const Float xyz[3], Float rgb[3]) {
    // Let R(lambda), B(lambda), and G(lambda) be spectral response curves and X(lambda),
    // Y(lambda), and Z(lambda) be the CIE XYZ color-matching curves. Let I be the integral
    // symbol.
    //
    // Then the input xyz vector is mapped to a unique rgb vector via matrix multiplication:
    // 
    // [r]   [I RX  I RY  I RZ][x]
    // [g] = [I GX  I GY  I GZ][y]
    // [b]   [I BX  I BY  I BZ][z]
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}

// Linearly maps RGB SPD coefficients to XYZ SPD coefficients.
inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
    // Let R(lambda), B(lambda), and G(lambda) be spectral response curves and X(lambda),
    // Y(lambda), and Z(lambda) be the CIE XYZ color-matching curves. Let I be the integral
    // symbol.
    //
    // Then the input rgb vector is mapped to a unique xyz vector via matrix multiplication:
    // 
    //                         -1
    // [x]   [I RX  I RY  I RZ]  [r]
    // [y] = [I GX  I GY  I GZ]  [g]
    // [z]   [I BX  I BY  I BZ]  [b]
    xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
    xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
    xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}

template <int nSpectrumSamples> class CoefficientSpectrum {
public:
    static const int nSamples = nSpectrumSamples;

    CoefficientSpectrum(Float v = 0.f) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] = v;
        }
    }

    CoefficientSpectrum &operator=(const CoefficientSpectrum &sp) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] = sp.c[i];
        }
        return *this;
    }

    Float &operator[](int i) {
        return c[i];
    }

    Float operator[](int i) const {
        return c[i];
    }

    // Spectrum addition is defined as the coefficient-wise sum of the 2 spectra.
    CoefficientSpectrum &operator+=(const CoefficientSpectrum &sp2) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] += sp2.c[i];
        }
        return *this;
    }

    CoefficientSpectrum operator+(const CoefficientSpectrum &sp2) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] += sp2.c[i];
        }
        return ret;
    }

    // Spectrum subtraction is defined as the coefficient-wise subtraction of the 2 spectra.
    CoefficientSpectrum operator-(const CoefficientSpectrum &sp2) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] -= sp2.c[i];
        }
        return ret;
    }

    // Spectrum product is defined as the coefficient-wise product of the 2 spectra.
    CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp2) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] *= sp2.c[i];
        }
        return *this;
    }

    CoefficientSpectrum operator*(const CoefficientSpectrum &sp2) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] *= sp2.c[i];
        }
        return ret;
    }

    // Spectrum division is defined as the coefficient-wise division of the 2 spectra.
    CoefficientSpectrum operator/(const CoefficientSpectrum &sp2) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] /= sp2.c[i];
        }
        return ret;
    }

    // Scalar multiplication multiplies each coefficient by the scalar.
    CoefficientSpectrum &operator*=(Float s) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] *= s;
        }
        return *this;
    }

    CoefficientSpectrum operator*(Float s) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] *= s;
        }
        return ret;
    }

    friend inline CoefficientSpectrum operator*(Float s, const CoefficientSpectrum &sp) {
        return sp * s;
    }

    // Scalar division divides each coefficient by the scalar.
    CoefficientSpectrum operator/(Float s) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] /= s;
        }
        Assert(!ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum &operator/=(Float s) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] /= s;
        }
        return *this;
    }

    bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            if (c[i] != sp.c[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }

    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = -c[i];
        }
        return ret;
    }

    // Spectrum square root is defined as the coefficient-wise square root of the spectrum.
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &sp) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = std::sqrt(sp.c[i]);
        }
        return ret;
    }

    // Spectrum power is defined as the coefficient-wise eth power of the spectrum.
    template <int n> friend
    inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n> &sp, Float e);

    // Spectrum exponentiation is defined as the coefficient-wise natural exponentiation of
    // the spectrum.
    friend CoefficientSpectrum Exp(const CoefficientSpectrum &sp) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = std::exp(sp.c[i]);
        }
        return ret;
    }

    CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = ::Clamp(c[i], low, high);
        }
        return ret;
    }

    bool IsBlack() const {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            if (c[i] != 0.) {
                return false;
            }
        }
        return true;
    }

    bool HasNaNs() const {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            if (std::isnan(c[i])) {
                return true;
            }
        }
        return false;
    }

    Float MaxComponentValue() const {
        Float m = c[0];
        for (int i = 1; i < nSpectrumSamples; ++i) {
            m = std::max(m, c[i]);
        }
        return m;
    }

protected:
    // Coefficients (of the linear combination of basis SPD functions that represents this SPD?).
    Float c[nSpectrumSamples];
};

// Spectrum power is defined as the coefficient-wise eth power of the spectrum.
template <int nSpectrumSamples>
inline CoefficientSpectrum<nSpectrumSamples> Pow(
    const CoefficientSpectrum<nSpectrumSamples> &sp, 
    Float e
) {
    CoefficientSpectrum<nSpectrumSamples> ret;
    for (int i = 0; i < nSpectrumSamples; ++i) {
        ret.c[i] = std::pow(s.c[i], e);
    }
    Assert(!ret.HasNaNs());
    return ret;
}

class RGBSpectrum : public CoefficientSpectrum<3> {
public:
    RGBSpectrum(Float v = 0.f) : CoefficientSpectrum<3>(v) {}

    RGBSpectrum(const CoefficientSpectrum<3> &v) : CoefficientSpectrum<3>(v) {}

    RGBSpectrum(const RGBSpectrum &s, SpectrumType type = SpectrumType::Reflectance) {
        *this = s;
    }

    static RGBSpectrum FromRGB(const Float rgb[3], SpectrumType = SpectrumType::Reflectance) {
        RGBSpectrum rsp;
        rsp.c[0] = rgb[0];
        rsp.c[1] = rgb[1];
        rsp.c[2] = rgb[2];
        return rsp;
    }

    void ToRGB(Float *rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }

    const RGBSpectrum &ToRGBSpectrum() const {
        return *this;
    }

    static RGBSpectrum FromXYZ(const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum rsp;
        XYZToRGB(xyz, rsp.c);
        return rsp;
    }

    void ToXYZ(Float xyz[3]) const {
        RGBToXYZ(c, xyz);
    }

    Float y() const {
        // Corresponds to the xyz[1] value returned by ToXYZ(xyz).
        const Float yWeight[3] = {0.212671f, 0.715160f, 0.072169f};
        return yWeight[0]*c[0] + yWeight[1]*c[1] + yWeight[2]*c[2];
    }

    static RGBSpectrum FromSampled(const Float *lambda, const Float *values, int n) {
        // Sort samples in order of increasing wavelength.
        if (!SpectrumSamplesSorted(lambda, n)) {
            std::vector<Float> sortedLambda(&lambda[0], &lambda[n]);
            std::vector<Float> sortedValues(&values[0], &values[n]);
            SortSpectrumSamples(&sortedLambda[0], &sortedValues[0], n);
            return FromSampled(&sortedLambda[0], &sortedValues[0], n);
        }

        // Map the sampled SPD to XYZ first and then to RGB.

        Float xyz[3] = {0.f, 0.f, 0.f};
        for (int i = 0; i < nCIESamples; ++i) {
            // Since the input SPD and the CIE XYZ color-matching curves may not be sampled
            // at the same wavelength intervals, we need to resample the SPD at the wavelength
            // intervals of the XYZ color-matching curves to match the samples with their
            // corresponding XYZ curve values.
            Float value = InterpolateSpectrumSamples(lambda, values, n, CIE_lambda[i]);

            // Sum of the products of this spectrum's SPD and the XYZ color-matching curves.
            xyz[0] += value * CIE_X[i];
            xyz[1] += value * CIE_Y[i];
            xyz[2] += value * CIE_Z[i];
        }

        Float lambdaDelta = Float(CIE_lambda[nCIESamples-1] - CIE_lambda[0]) / Float(nCIESamples);

        // Riemman sums that approximate the definite integrals of the products of the SPD
        // and the XYZ color-matching curves.
        xyz[0] *= lambdaDelta;
        xyz[1] *= lambdaDelta;
        xyz[2] *= lambdaDelta;

        // ?
        xyz[0] *= 1.f / CIE_Y_integral;
        xyz[1] *= 1.f / CIE_Y_integral;
        xyz[2] *= 1.f / CIE_Y_integral;

        return FromXYZ(xyz);
    }
};

class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
private:
    // n samples of the CIE X(lambda), Y(lambda), and Z(lambda) 
    // color-matching curves, where n = nSpectralSamples.
    //
    // The sampled XYZ curves CIE_X, CIE_Y, and CIE_Z are sampled at 1nm intervals.
    // This SampledSpectrum is sampled at (sampledLambdaEnd-sampledLambdaStart)/nSpectralSamples
    // nm intervals, for which the CIE_X, CIE_Y, and CIE_Z curve arrays may not have samples. So
    // the average of the CIE_X, CIE_Y, and CIE_Z samples contained by one of this SampledSpectrum's
    // intervals is computed and regarded as the CIE X, Y, Z value of that interval. 
    static SampledSpectrum X;
    static SampledSpectrum Y;
    static SampledSpectrum Z;

    // SPDs for RGB primaries and colors that are mixes of primaries exclusively.
    // All these SPDs exist at the namespace scope with 32 samples. This SampledSpectrum resamples
    // the same SPDs at n = nSpectralSamples samples.
    static SampledSpectrum rgbRefl2SpectWhite;
    static SampledSpectrum rgbRefl2SpectCyan;
    static SampledSpectrum rgbRefl2SpectMagenta;
    static SampledSpectrum rgbRefl2SpectYellow;
    static SampledSpectrum rgbRefl2SpectRed;
    static SampledSpectrum rgbRefl2SpectGreen;
    static SampledSpectrum rgbRefl2SpectBlue;
    static SampledSpectrum rgbIllum2SpectWhite;
    static SampledSpectrum rgbIllum2SpectCyan;
    static SampledSpectrum rgbIllum2SpectMagenta;
    static SampledSpectrum rgbIllum2SpectYellow;
    static SampledSpectrum rgbIllum2SpectRed;
    static SampledSpectrum rgbIllum2SpectGreen;
    static SampledSpectrum rgbIllum2SpectBlue;
    
public:
    SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) {}

    SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v) 
        : CoefficientSpectrum<nSpectralSamples>(v) {}

    SampledSpectrum(const RGBSpectrum &rsp, SpectrumType type) {
        Float rgb[3];
        rsp.ToRGB(rgb);
        *this = SampledSpectrum::FromRGB(rgb, type);
    };

    static void Init() {
        // Sample the XYZ color-matching curves at (sampledLambdaEnd-sampledLambdaStart)/nSpectralSamples
        // nm intervals.
        for (int i = 0; i < nSpectralSamples; ++i) {
            Float lambda0 = Lerp(
                Float(i) / Float(nSpectralSamples),
                sampledLambdaStart, 
                sampledLambdaEnd
            );

            Float lambda1 = Lerp(
                Float(i+1) / Float(nSpectralSamples),
                sampledLambdaStart,
                sampledLambdaEnd
            );

            // The XYZ value of the interval is the average of the CIE_X, CIE_Y, and CIE_Z samples
            // contained by the interval.
            X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, lambda0, lambda1);
            Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, lambda0, lambda1);
            Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, lambda0, lambda1);
        }

        // Resample the (namespace-level) precomputed RGB SPDs at n = nSpectralSamples samples.
        for (int i = 0; i < nSpectralSamples; ++i) {
            Float lambda0 = Lerp(Float(i) / nSpectralSamples, sampledLambdaStart, sampledLambdaEnd);
            Float lambda1 = Lerp(Float(i+1) / nSpectralSamples, sampledLambdaStart, sampledLambdaEnd);

            rgbRefl2SpectWhite.c[i]   = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectCyan.c[i]    = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectYellow.c[i]  = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectRed.c[i]     = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectGreen.c[i]   = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen, nRGB2SpectSamples, lambda0, lambda1);
            rgbRefl2SpectBlue.c[i]    = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue, nRGB2SpectSamples, lambda0, lambda1);

            rgbIllum2SpectWhite.c[i]   = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectCyan.c[i]    = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectYellow.c[i]  = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectRed.c[i]     = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectGreen.c[i]   = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen, nRGB2SpectSamples, lambda0, lambda1);
            rgbIllum2SpectBlue.c[i]    = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue, nRGB2SpectSamples, lambda0, lambda1);
        }
    }

    // n (wavelength, value) pairs.
    static SampledSpectrum FromSampled(const Float *lambda, const Float *values, int n) {
        // Sort samples in order of increasing wavelength.
        if (!SpectrumSamplesSorted(lambda, n)) {
            std::vector<Float> sortedLambda(&lambda[0], &lambda[n]);
            std::vector<Float> sortedValues(&values[0], &values[n]);
            SortSpectrumSamples(&sortedLambda[0], &sortedValues[0], n);
            return FromSampled(&sortedLambda[0], &sortedValues[0], n);
        }

        SampledSpectrum r;
        for (int i = 0; i < nSpectralSamples; ++i) {
            // i and i+1 are the endpoints of the ith wavelength subinterval.
            // lambda0 and lambda1 are the wavelengths of these endpoints.
            Float lambda0 = Lerp(Float(i) / nSpectralSamples, Float(sampledLambdaStart), Float(sampledLambdaEnd));
            Float lambda1 = Lerp(Float(i+1) / nSpectralSamples, Float(sampledLambdaStart), Float(sampledLambdaEnd));

            r.c[i] = AverageSpectrumSamples(lambda, values, n, lambda0, lambda1);
        }
        return r;
    }

    // Maps this spectrum's SPD to its equivalent XYZ tristimulus value.
    //
    // X is defined as the definite integral of the product of the SPD and the CIE X curve,
    // and is approximated by the corresponding Riemann sum over n=nSpectralSamples subintervals.
    //
    // Y and Z are defined analogously. 
    void ToXYZ(Float xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = 0.f;

        // Sum of the products of this spectrum's SPD and the XYZ color-matching curves.
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += c[i] * X.c[i];
            xyz[1] += c[i] * Y.c[i];
            xyz[2] += c[i] * Z.c[i];
        }

        Float lambdaDelta = Float(sampledLambdaEnd-sampledLambdaStart) / Float(nSpectralSamples);

        // Riemman sums that approximate the definite integrals of the products of the SPD 
        // and the XYZ color-matching curves.
        xyz[0] *= lambdaDelta;
        xyz[1] *= lambdaDelta;
        xyz[2] *= lambdaDelta;

        // ?
        xyz[0] *= 1.f / CIE_Y_integral;
        xyz[1] *= 1.f / CIE_Y_integral;
        xyz[2] *= 1.f / CIE_Y_integral;
    }

    // XYZ y tristimulus value.
    Float y() const {
        Float yy = 0.f;

        for (int i = 0; i < nSpectralSamples; ++i) {
            yy += c[i] * Y.c[i];
        }

        Float lambdaDelta = Float(sampledLambdaEnd-sampledLambdaStart) / Float (nSpectralSamples);

        // Riemann sum.
        yy *= lambdaDelta;

        // ?
        return yy / CIE_Y_integral;
    }

    // Maps this SPD's coefficients to RGB coefficients by first mapping them to XYZ.
    void ToRGB(Float rgb[3]) const {
        Float xyz[3];

        // This spectrum's XYZ coefficients.
        ToXYZ(xyz);

        // Map this spectrum's XYZ coefficients to RGB.
        XYZToRGB(xyz, rgb);
    }

    RGBSpectrum ToRGBSpectrum() const {
        Float rgb[3];
        ToRGB(rgb);
        return RGBSpectrum::FromRGB(rgb);
    }

    static SampledSpectrum FromRGB(const Float rgb[3], SpectrumType type);

    static SampledSpectrum FromXYZ(const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
        Float rgb[3];
        XYZToRGB(xyz, rgb);
        return FromRGB(rgb, type);
    }
};

inline Spectrum Lerp(Float t, const Spectrum &sp1, const Spectrum &sp2) {
    return (1 - t)*sp1 + t*sp2;
}

#endif // CPBRT_CORE_SPECTRUM_H