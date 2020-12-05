#include <cmath>

#include "cpbrt.h"

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


inline Spectrum Lerp(Float t, const Spectrum &sp1, const Spectrum &sp2) {
    return (1 - t)*sp1 + t*sp2;
}