#include "mipmap.h"

template <typename T> const T &MIPMap<T>::Texel(int level, int s, int t) const {
    const BlockedArray<T> &l = *pyramid[level];

    switch (wrapMode) {
        case ImageWrap::Repeat:
            s = Mod(s, l.uSize());
            t = Mod(t, l.vSize());
            break;
        case ImageWrap::Clamp:
            s = Clamp(s, 0, l.uSize() - 1);
            t = Clamp(t, 0, l.vSize() - 1);
            break;
        case ImageWrap::Black: {
            static const T black = 0.f;
            if (s < 0 || s >= (int)l.uSize() || t < 0 || t >= (int)l.vSize()) {
                return black;
            }
            break;
        }
    }

    return l(s, t);
}

template <typename T> T MIPMap<T>::Lookup(const Point2f &st, Float width) const {
    // Compute mipmap level. The level may be a floating-point number between 2
    // pyramid levels. If i < level < i+1 are levels, then (i - level) measures how
    // close level is to pyramid levels i and i+1. (i - level) serves then as interpolator
    // between the filtered sample from level i and the filtered sample from level i+1.
    Float level = Levels() - 1 + Log2(std::max(width, (Float) 1e-8));

    // Perform trilinear interpolation at computed mipmap level.
    if (level < 0) {
        return triangle(0, st);
    } else if (level >= Levels() - 1) {
        return Texel(Levels() - 1, 0, 0);
    } else {
        int iLevel = std::floor(level);
        Float delta = level - iLevel;
        return Lerp(delta, triangle(iLevel, st), triangle(iLevel + 1, st));
    }
}