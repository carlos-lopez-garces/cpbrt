#ifndef CPBRT_CORE_EFLOAT_H
#define CPBRT_CORE_EFLOAT_H

#include "cpbrt.h"
#include "stringprint.h"

// Represents the result of a floating-point arithmetic operation as an interval
// [low, high] of possible values, among which is the actual computed floating-point number.
// For that computed floating-point number, the real number may be any of [low, high].
// This is called interval arithmetic and is a rounding error analysis technique.
class EFloat {
public:
    EFloat() {}

    EFloat(float v, float err = 0.f) : v(v) {
        if (err == 0.) {
            low = high = v;
        } else {
            // Lower error bound.
            low = NextFloatDown(v - err);
            // Upper error bound.
            high = NextFloatUp(v + err);
        }
    }

    EFloat(const EFloat &ef) {
        v = ef.v;
        low = ef.low;
        high = ef.high;
    }

    float UpperBound() const {
        return high;
    }

    float LowerBound() const {
        return low;
    }

    // Explicit cast to float.
    explicit operator float() const {
        return v;
    }

    // Explicit cast to double.
    explicit operator double() const {
        return v;
    }

     EFloat &operator=(const EFloat &ef) {
        if (&ef != this) {
            v = ef.v;
            low = ef.low;
            high = ef.high;
        }

        return *this;
    }

    inline bool operator==(EFloat fe) const {
        return v == fe.v;
    }

    EFloat operator+(EFloat ef) const {
        EFloat r;
        
        r.v = v + ef.v;

        // Running sum of lower error bound.
        r.low = NextFloatDown(LowerBound() + ef.LowerBound());

        // Running sum of upper error bound.
        r.high = NextFloatUp(UpperBound() + ef.UpperBound());

        return r;
    }

    EFloat operator-(EFloat ef) const {
        EFloat r;
        
        r.v = v - ef.v;

        // While sum adds to the error, subtraction reduces it.

        // Running sum of lower error bound.
        r.low = NextFloatDown(LowerBound() - ef.UpperBound());

        // Running sum of upper error bound.
        r.high = NextFloatUp(UpperBound() - ef.LowerBound());

        return r;
    }

    EFloat operator*(EFloat ef) const {
        EFloat r;
        
        r.v = v * ef.v;

        Float prod[4] = {
            LowerBound() * ef.LowerBound(),
            UpperBound() * ef.LowerBound(),
            LowerBound() * ef.UpperBound(),
            UpperBound() * ef.UpperBound()
        };

        // TODO: explain.
        r.low = NextFloatDown(
            std::min(
                std::min(prod[0], prod[1]),
                std::min(prod[2], prod[3])
            )
        );

        r.high = NextFloatUp(
            std::max(
                std::max(prod[0], prod[1]),
                std::max(prod[2], prod[3])
            )
        );

        return r;
    }

    EFloat operator/(EFloat ef) const {
        EFloat r;

        r.v = v / ef.v;

        // TODO: explain.
        if (ef.low < 0 && ef.high > 0) {
            r.low = -Infinity;
            r.high = Infinity;
        } else {
            Float div[4] = {
                LowerBound() / ef.LowerBound(),
                UpperBound() / ef.LowerBound(),
                LowerBound() / ef.UpperBound(),
                UpperBound() / ef.UpperBound()
            };

            r.low = NextFloatDown(
                std::min(
                    std::min(div[0], div[1]),
                    std::min(div[2], div[3])
                )
            );

            r.high = NextFloatUp(
                std::max(
                    std::max(div[0], div[1]),
                    std::max(div[2], div[3])
                )
            );
        }

        return r;
    }

    EFloat operator-() const {
        EFloat r;

        r.v = -v;

        r.low = -high;
        r.high = -low;
        
        return r;
    }

    friend std::ostream &operator<<(std::ostream &os, const EFloat &ef) {
        os << StringPrintf("v=%f (%a) - [%f, %f]", ef.v, ef.v, ef.low, ef.high);
        return os;
    }

    float GetAbsoluteError() const {
        // Always be conservative. Return the largest of the 2 bounds and round up.
        return NextFloatUp(std::max(std::abs(high - v), std::abs(v - low)));
    }

    // Friend functions.
    friend inline EFloat sqrt(EFloat fe);
    friend inline EFloat abs(EFloat fe);
    // Solves the quadratic equation represented by the given coefficients: At^2 + Bt + C = 0.
    // If 2 solutions exist, the smallest one is returned in t0.
    friend inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1);

private:
    // Computed value.
    float v;
    float low;
    float high;
};

// The other operand being a regular float.

inline EFloat operator+(float f, EFloat fe) {
    return EFloat(f) + fe;
}

inline EFloat operator/(float f, EFloat fe) {
    return EFloat(f) / fe; 
}

inline EFloat operator*(float f, EFloat fe) {
    return EFloat(f) * fe;
}

inline EFloat operator-(float f, EFloat fe) { 
    return EFloat(f) - fe; 
}

inline EFloat sqrt(EFloat fe) {
    EFloat r;

    r.v = std::sqrt(fe.v);

    // Square root reduces the error.
    r.low = NextFloatDown(std::sqrt(fe.low));
    r.high = NextFloatUp(std::sqrt(fe.high));

    return r;
}

inline EFloat abs(EFloat fe) {
    if (fe.low >= 0) {
        // The entire interval is greater than zero; they are already their absolute values.
        return fe;
    }
    else if (fe.high <= 0) {
        // The entire interval is less than zero.
        EFloat r;

        r.v = -fe.v;

        r.low = -fe.high;
        r.high = -fe.low;

        return r;
    } else {
        // The interval straddles zero.
        EFloat r;

        r.v = std::abs(fe.v);

        r.low = 0;
        r.high = std::max(-fe.low, fe.high);

        return r;
    }
}

inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1) {
    // Discriminant.
    double discrim = (double) B.v * (double) B.v - 4. * (double) A.v * (double) C.v;
    if (discrim < 0.) {
        return false;
    }
    double rootDiscrim = std::sqrt(discrim);

    EFloat floatRootDiscrim(rootDiscrim, MachineEpsilon * rootDiscrim);

    // Compute quadratic _t_ values
    EFloat q;
    if ((float)B < 0) {
        q = -.5 * (B - floatRootDiscrim);
    } else {
        q = -.5 * (B + floatRootDiscrim);
    }

    *t0 = q / A;
    *t1 = C / q;
    if ((float)*t0 > (float)*t1) {
        std::swap(*t0, *t1);
    }

    return true;
}

#endif // CPBRT_CORE_EFLOAT_H