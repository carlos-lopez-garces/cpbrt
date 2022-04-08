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