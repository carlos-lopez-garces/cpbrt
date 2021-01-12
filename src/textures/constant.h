#include "core/interaction.h"
#include "core/texture.h"

template <typename T> class ConstantTexture : public Texture<T> {
private:
    // Constant.  
    T value;

public:
    ConstantTexture(const T &value) : value(value) {}

    T Evaluate(const SurfaceInteraction &) const {
        return value;
    }
};