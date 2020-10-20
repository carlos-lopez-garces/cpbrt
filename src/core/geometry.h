#include "pbrt.h"

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
template <typename T> inline Vector2<T> Max(const Vector2<T> &v1, const Vector2<T> &v2) {
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

    // Cast Normal3<T> to Vector3<T>.
    explicit Vector3(const Normal3<T> &n);

    Vector3<T> operator+(const Vector3<T> &v) const {
        return Vector3<T>(x+v.x, y+v.y, z+v.z);
    }

    Vector3<T> operator-(const Vector3<T> &v) const {
        return Vector3<T>(x-v.x, y-v.y, z-v.z);
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
template <typename T> inline Vector3<T> Max(const Vector3<T> &v1, const Vector3<T> &v2) {
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

template <typename T> class Point2 {
public:
    T x, y;

    Point2() {
        x = y = 0;
    }

    Point2(T x, T y) : x(x), y(y) {
        Assert(!HasNaNs());
    }

    // Require explicit cast to convert from Point3 to Point2.
    explicit Point2(const Point3<T> &p) : x(p.x), y(p.y) {
        Assert(!HasNaNs());
    }

    // Require explicit type to convert from Point2<U> to Point2<T>.
    template <typename U> explicit Point2(const Point2<U> &p)
        : x((T)p.x), y((T)p.y)
    {
        Assert(!HasNaNs());
    }

    // Cast Point2<T> to Vector2<U>.
    template <typename U> explicit operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }

    Point2<T> operator+(const Vector2<T> &v) const {
        return Point2<T>(x+v.x, y+v.y);
    }

    Point2<T> &operator+=(const Vector2<T> &v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    // Addition of points only makes sense in the context of computing a weighted
    // sum of points, with the sum of the weights being 1.
    Point2<T> operator+(const Point2<T> &p) const {
        return Point3<T>(x+p.x, y+p.y);
    }

    // Subtracting a point from this point results in the vector between them.
    Vector2<T> operator-(const Point2<T> &p) const {
        return Vector2<T>(x-p.x, y-p.y);
    }

    // Subtracting a vector from this point results in a displaced point.
    Point2<T> operator-(const Vector2<T> &v) const {
        return Point2<T>(x-v.x, y-v.y);
    }

    // Subtracting a vector from this point translates this point.
    Point2<T> &operator-=(const Vector2<T> &v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    // Scalar multiplication of points only makes sense in the context of 
    // computing a weighted sum of points, with the sum of the weights being 1.
    Point2<T> operator*(T s) const {
        return Point2<T>(s*x, s*y);
    }

    bool HasNaNs() {
        return std::isnan(x) || std::isnan(y);
    }
};

typedef Point2<Float> Point2f;
typedef Point2<int> Point2i;

template <typename T> inline Float Distance(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1-p2).Length();
}

template <typename T> inline Float DistanceSquared(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1-p2).LengthSquared();
}

template <typename T> inline Point2<T> operator*(T s, const Point2<T> &p) {
    return p * s;
}

// Linear interpolation (0 <= t <= 1) and extrapolation (t < 0, t > 1).
template <typename T> inline Point2<T> Lerp(Float t, const Point2<T> &p0, const Point2<T> &p1) {
    return (1-t)*p0 + t*p1;
}

// Component-wise minimum.
template <typename T> inline Point2<T> Min(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(
        std::min(p1.x, p2.x),
        std::min(p1.y, p2.y),
    );
}

// Component-wise maximum.
template <typename T> inline Point2<T> Max(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(
        std::max(p1.x, p2.x),
        std::max(p1.y, p2.y),
    );
}

template <typename T> inline Point2<T> Floor(const Point2<T> &p) {
    return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T> inline Point2<T> Ceil(const Point2<T> &p) {
    return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T> inline Point2<T> Abs(const Point2<T> &p) {
    return Point2<T>(std::abs(p.x), std::abs(p.y));
}

template <typename T> inline Point2<T> Permute(const Point2<T> &p, int x, int y) {
    return Point2<T>(p[x], p[y]);
}

template <typename T> class Point3 {
public:
    T x, y, z;

    Point3() {
        x = y = z = 0;
    }

    Point3(T x, T y, T z) : x(x), y(y), z(z) {
        Assert(!HasNaNs());
    }

    // Require explicit type to convert from Point3<U> to Point3<T>.
    template <typename U> explicit Point3(const Point3<U> &p) 
        : x((T)p.x), y((T)p.y), z((T)p.z) 
    {
        Assert(!HasNaNs());
    }

    // Cast Point3<T> to Vector3<U>.
    template <typename U> explicit operator Vector3<U>() const {
        return Vector3<U>(x, y, z);
    }

    Point3<T> operator+(const Vector3<T> &v) const {
        return Point3<T>(x+v.x, y+v.y, z+v.z);
    }

    Point3<T> &operator+=(const Vector3<T> &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    // Addition of points only makes sense in the context of computing a weighted
    // sum of points, with the sum of the weights being 1.
    Point3<T> operator+(const Point3<T> &p) const {
        return Point3<T>(x+p.x, y+p.y, z+p.z);
    }

    // Subtracting a point from this point results in the vector between them.
    Vector3<T> operator-(const Point3<T> &p) const {
        return Vector3<T>(x-p.x, y-p.y, z-p.z);
    }

    // Subtracting a vector from this point results in a displaced point.
    Point3<T> operator-(const Vector3<T> &v) const {
        return Point3<T>(x-v.x, y-v.y, z-v.z);
    }

    // Subtracting a vector from this point translates this point.
    Point3<T> &operator-=(const Vector3<T> &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    // Scalar multiplication of points only makes sense in the context of computing a weighted
    // sum of points, with the sum of the weights being 1.
    Point3<T> operator*(T s) const {
        return Point3<T>(s*x, s*y, s*z);
    }

    bool HasNaNs() {
        return std::isnan(x) || std::isnan(y);
    }
};

typedef Point3<Float> Point3f;
typedef Point3<int> Point3i;

template <typename T> inline Float Distance(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1-p2).Length();
}

template <typename T> inline Float DistanceSquared(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1-p2).LengthSquared();
}

template <typename T> inline Point3<T> operator*(T s, const Point3<T> &p) {
    return p * s;
}

// Linear interpolation (0 <= t <= 1) and extrapolation (t < 0, t > 1).
template <typename T> inline Point3<T> Lerp(Float t, const Point3<T> &p0, const Point3<T> &p1) {
    return (1-t)*p0 + t*p1;
}

// Component-wise minimum.
template <typename T> inline Point3<T> Min(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(
        std::min(p1.x, p2.x),
        std::min(p1.y, p2.y),
        std::min(p1.z, p2.z)
    );
}

// Component-wise maximum.
template <typename T> inline Point3<T> Max(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(
        std::max(p1.x, p2.x),
        std::max(p1.y, p2.y),
        std::max(p1.z, p2.z)
    );
}

template <typename T> inline Point3<T> Floor(const Point3<T> &p) {
    return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T> inline Point3<T> Ceil(const Point3<T> &p) {
    return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T> inline Point3<T> Abs(const Point3<T> &p) {
    return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T> inline Point3<T> Permute(const Point3<T> &p, int x, int y, int z) {
    return Point3<T>(p[x], p[y], p[z]);
}

// Normals are not necessarily normalized.
template <typename T> class Normal3 {
public:
    T x, y, z;

    Normal3() {
        x = y = z = 0;
    }

    Normal3(T x, T y, T z) : x(x), y(y), z(z) {
        Assert(!HasNans());
    }

    explicit Normal3<T>(const Vector3<T> &v)
        : x(v.x), y(v.y), z(v.z) 
    {
        Assert(!HasNaNs());
    }

    Normal3<T> operator+(const Normal3<T> &n) const {
        return Normal3<T>(x+n.x, y+n.y, z+n.z);
    }

    Normal3<T> operator-(const Normal3<T> &n) const {
        return Normal3<T>(x-n.x, y-n.y, z-n.z);
    }

    Normal3<T> operator*(T s) const {
        return Normal3<T>(s*x, s*y, s*z);
    }

    Normal3<T> operator-() const {
        return Normal3<T>(-x, -y, -z);
    }

    Float LengthSquared() const {
        return x*x + y*y + z*z;
    }

    Float Length() const {
        return std::sqrt(LengthSquared());
    }

    bool HasNaNs() {
        return std::isnan(x) || std::isnan(y);
    }
};

typedef Normal3<Float> Normal3f;

template <typename T> inline Normal3<T> operator*(T s, const Normal3<T> &n) {
    return n * s;
}

// Normals are not necessarily normalized.
template <typename T> inline Normal3<T> Normalize(const Normal3<T> &n) {
    return n / n.Length();
}

template <typename T> inline Vector3<T>::Vector3(const Normal3<T> &n)
    : x(n.x), y(n.y), z(n.z)
{
    Assert(!HasNaNs());
}

template <typename T> inline T Dot(const Normal3<T> &n, const Vector3<T> &v) {
    return n.x*v.x + n.y*v.y + n.z*v.z;
}

template <typename T> inline T AbsDot(const Normal3<T> &n, const Vector3<T> &v) {
    return std::abs(Dot(n, v));
}

template <typename T> inline T Dot(const Vector3<T> &v, const Normal3<T> &n) {
    return Dot(v, n);
}

template <typename T> inline T AbsDot(const Vector3<T> &v, const Normal3<T> &n) {
    return std::abs(Dot(v, n));
}

// Flip a normal so that it lies in the same hemisphere as the vector.
template <typename T> inline Normal3<T> FaceForward(
    const Normal3<T> &n, const Vector3<T> &v
) {
    return (Dot(n, v) < 0.f) ? -n : n;
}