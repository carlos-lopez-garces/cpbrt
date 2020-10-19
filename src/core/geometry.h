template <typename T> class Vector2 {
public:
    T x, y;

    Vector2() {
        x = y = 0;
    }

    Vector2(T x, T y) : x(x), y(y) {
        Assert(!HasNaNs());
    }

    Vector2<T> operator+(const Vector2<T> &v) const {
        return Vector2<T>(x+v.x, y+v.y);
    }

    Vector2<T> operator-(const Vector2<T> &v) const {
        return Vector2<T>(x-v.x, y-v.y);
    }

    Vector2<T> operator*(T s) const {
        return Vector2<T>(s*x, s*y);
    }

    Vector2<T> operator/(T f) const {
        Assert(f != 0);
        // Multiplication is faster than division. Compilers are generally restricted
        // from performing transformations of this type. So they wouldn't happen automatically.
        Float reciprocal = (Float)1/f;
        return Vector2<T>(x*reciprocal, y*reciprocal);
    }

    Vector2<T> operator-() const {
        return Vector2<T>(-x, -y);
    }

    Vector2<T> &operator+=(const Vector2<T> &v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    Vector2<T> &operator-=(const Vector2<T> &v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    Vector2<T> &operator*=(T s) {
        x *= s;
        y *= s;
        return *this;
    }

    Vector2<T> &operator/=(T f) {
        Assert(f != 0);
        Float reciprocal = (Float)1/f;
        x *= reciprocal;
        y *= reciprocal;
        return *this;
    }

    Float LengthSquared() const {
        return x*x + y*y;
    }

    Float Length() const {
        return std::sqrt(LengthSquared());
    }

    T operator[](int i) const {
        Assert(i >= 0 && i <= 1);
        return (i == 0) ? x : y;
    }

    T &operator[](int i) const {
        Assert(i >= 0 && i <= 1);
        return (i == 0) ? x : y;
    }

    bool HasNaNs() {
        return std::isnan(x) || std::isnan(y);
    }
};

typedef Vector2<Float> Vector2f;
typedef Vector2<int> Vector2i;

template <typename T> inline Vector2<T> operator*(T s, const Vector2<T> &v) {
    return v * s;
}

template <typename T> inline Vector2<T> Abs(const Vector2<T> &v) {
    return Vector2<T>(std::abs(v.x), std::abs(v.y));
}

template <typename T> inline T Dot(const Vector2<T> &v1, const Vector2<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y;
}

template <typename T> inline T AbsDot(const Vector2<T> &v1, const Vector2<T> &v2) {
    return std::abs(Dot(v1, v2));
}

template <typename T> inline Vector2<T> Normalize(const Vector2<T> &v) {
    return v / v.Length();
}

template <typename T> inline T MinComponent(const Vector2<T> &v) {
    return std::min(v.x, v.y);
}

template <typename T> inline T MaxComponent(const Vector2<T> &v) {
    return std::max(v.x, v.y);
}

template <typename T> inline int MaxDimension(const Vector2<T> &v) {
    return (v.x > v.y) ? 0 : 1;
}

// Component-wise minimum.
template <typename T> inline Vector2<T> Min(const Vector2<T> &v1, const Vector2<T> &v2) {
    return Vector2<T>(
        std::min(v1.x, v2.x),
        std::min(v1.y, v2.y),
    );
}

// Component-wise maximum.
template <typename T> inline Vector2<T> Min(const Vector2<T> &v1, const Vector2<T> &v2) {
    return Vector2<T>(
        std::max(v1.x, v2.x),
        std::max(v1.y, v2.y),
    );
}

// Permutes the coordinate values according to the index values provided.
template <typename T> inline Vector2<T> Permute(const Vector2<T> &v, int x, int y) {
    return Vector2<T>(v[x], v[y]);
}

template <typename T> class Vector3 {
public:
    T x, y, z;

    Vector3() {
        x = y = z = 0;
    }

    Vector3(T x, T y, T z) : x(x), y(y), z(z) {
        Assert(!HasNans());
    }

    Vector3<T> operator+(const Vector3<T> &v) const {
        return Vector3<T>(x+v.x, y+x.y, z+v.z);
    }

    Vector3<T> operator-(const Vector3<T> &v) const {
        return Vector3<T>(x-v.x, y-x.y, z-v.z);
    }

    Vector3<T> operator*(T s) const {
        return Vector3<T>(s*x, s*y, s*z);
    }

    Vector3<T> operator/(T f) const {
        Assert(f != 0);
        // Multiplication is faster than division. Compilers are generally restricted
        // from performing transformations of this type. So they wouldn't happen automatically.
        Float reciprocal = (Float)1/f;
        return Vector3<T>(x*reciprocal, y*reciprocal, z*reciprocal);
    }

    Vector3<T> operator-() const {
        return Vector3<T>(-x, -y, -z);
    }

    Vector<T> &operator+=(const Vector3<T> &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vector<T> &operator-=(const Vector3<T> &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    Vector3<T> &operator*=(T s) {
        x *= s;
        y *= s;
        z *= z;
        return *this;
    }

    Vector3<T> &operator/=(T f) {
        Assert(f != 0);
        Float reciprocal = (Float)1/f;
        x *= reciprocal;
        y *= reciprocal;
        z *= reciprocal;
        return *this;
    }

    Float LengthSquared() const {
        return x*x + y*y + z*z;
    }

    Float Length() const {
        return std::sqrt(LengthSquared());
    }

    T operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (i == 0) 
            ? x 
            : (i == 1) ? y : z;
    }

    T &operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        return (i == 0) 
            ? x 
            : (i == 1) ? y : z;
    }

    bool HasNaNs() {
        return std::isnan(x) || std::isnan(y) || std::isnan(z);
    }
};

typedef Vector3<Float> Vector3f;
typedef Vector3<int> Vector3i;

template <typename T> inline Vector3<T> operator*(T s, const Vector3<T> &v) {
    return v * s;
}

template <typename T> inline Vector3<T> Abs(const Vector3<T> &v) {
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T> inline T Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template <typename T> inline T AbsDot(const Vector3<T> &v1, const Vector3<T> &v2) {
    return std::abs(Dot(v1, v2));
}

template <typename T> inline Vector3<T> Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
    // Use double precision to avoid _catastrophic cancellation_.
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>(
        (v1y*v2z) - (v1z*v2y),
        (v1z*v2x) - (v1x*v2z),
        (v1x*v2y) - (v1y*v2x)
    );
}

template <typename T> inline Vector3<T> Normalize(const Vector3<T> &v) {
    return v / v.Length();
}

template <typename T> inline T MinComponent(const Vector3<T> &v) {
    return std::min(v.x, std::min(v.y, v.z));
}

template <typename T> inline T MaxComponent(const Vector3<T> &v) {
    return std::max(v.x, std::max(v.y, v.z));
}

template <typename T> inline int MaxDimension(const Vector3<T> &v) {
    return (v.x > v.y) 
        ? ((v.x > v.z) ? 0 : 2) 
        : ((v.y > v.z) ? 1 : 2);
}

// Component-wise minimum.
template <typename T> inline Vector3<T> Min(const Vector3<T> &v1, const Vector3<T> &v2) {
    return Vector3<T>(
        std::min(v1.x, v2.x),
        std::min(v1.y, v2.y),
        std::min(v1.z, v2.z)
    );
}

// Component-wise maximum.
template <typename T> inline Vector3<T> Min(const Vector3<T> &v1, const Vector3<T> &v2) {
    return Vector3<T>(
        std::max(v1.x, v2.x),
        std::max(v1.y, v2.y),
        std::max(v1.z, v2.z)
    );
}

// Permutes the coordinate values according to the index values provided.
template <typename T> inline Vector3<T> Permute(const Vector3<T> &v, int x, int y, int z) {
    return Vector3<T>(v[x], v[y], v[z]);
}

// Assumes that v1 is normalized.
template <typename T> inline void CoordinateSystem(
    const Vector3<T> &v1, Vector3<T> *v2, Vector3<T> *v3
) {
    if (std::abs(v1.x) > std::abs(v1.y)) {
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x*v1.x + v1.z*v1.z);
    } else {
        *v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y*v1.y + v1.z*v1.z);
    }

    *v3 = Cross(v1, *v2);
}