#include <cmath>

#include "cpbrt.h"

static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;

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

inline Spectrum Lerp(Float t, const Spectrum &sp1, const Spectrum &sp2) {
    return (1 - t)*sp1 + t*sp2;
}

template <int nSpectrumSamples> class CoefficientSpectrum {
public:
    static const int nSamples = nSpectrumSamples;

    CoefficientSpectrum(Float v = 0.f) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            c[i] = v;
        }
    }

    Float &operator[](int i) {
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
        for (int = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] /= s;
        }
        Assert(!ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum &operator/=(Float s) {
        for (int = 0; i < nSpectrumSamples; ++i) {
            c[i] /= s;
        }
        return *this;
    }

    bool operator==(const CoefficientSpectrum &sp) const {
        for (int = 0; i < nSpectrumSamples; ++i) {
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
        for (int = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = -c[i];
        }
        return ret;
    }

    // Spectrum square root is defined as the coefficient-wise square root of the spectrum.
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &sp) {
        CoefficientSpectrum ret;
        for (int = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = std::sqrt(s.c[i]);
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
        for (int = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = std::exp(s.c[i]);
        }
        return ret;
    }

    CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const {
        CoefficientSpectrum ret;
        for (int = 0; i < nSpectrumSamples; ++i) {
            ret.c[i] = ::Clamp(c[i], low, high);
        }
        return ret;
    }

    bool IsBlack() const {
        for (int = 0; i < nSpectrumSamples; ++i) {
            if (c[i] != 0.) {
                return false;
            }
        }
        return true;
    }

    bool HasNaNs() const {
        for (int = 0; i < nSpectrumSamples; ++i) {
            if (std::isnan(c[i])) {
                return true;
            }
        }
        return false;
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
    for (int = 0; i < nSpectrumSamples; ++i) {
        ret.c[i] = std::pow(s.c[i], e);
    }
    Assert(!ret.HasNaNs());
    return ret;
}

class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
public:
    SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) {}

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
private:
};