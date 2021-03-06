#ifndef CPBRT_CORE_GEOMETRY_H
#define CPBRT_CORE_GEOMETRY_H

#include <cmath>
#include <iterator>

#include "cpbrt.h"

template <typename T> inline bool isNaN(const T x) {
    return std::isnan(x);
}

template <> inline bool isNaN(const int x) {
    return false;
}

template <typename T> class Vector2 {
public:
    T x, y;

    Vector2() {
        x = y = 0;
    }

    Vector2(T x, T y) : x(x), y(y) {
        Assert(!HasNaNs());
    }

    bool operator==(const Vector2<T> &v) const {
        return x == v.x && y == v.y;
    }

    bool operator!=(const Vector2<T> &v) const {
        return x != v.x || y != v.y;
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

    T &operator[](int i) {
        Assert(i >= 0 && i <= 1);
        return (i == 0) ? x : y;
    }

    bool HasNaNs() {
        return isNaN(x) || isNaN(y);
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
        Assert(!HasNaNs());
    }

    // Cast Normal3<T> to Vector3<T>.
    explicit Vector3(const Normal3<T> &n);

    bool operator==(const Vector3<T> &v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const Vector3<T> &v) const {
        return x != v.x || y != v.y || z != v.z;
    }

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

    Vector3<T> &operator+=(const Vector3<T> &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vector3<T> &operator-=(const Vector3<T> &v) {
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

    T &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        return (i == 0) 
            ? x 
            : (i == 1) ? y : z;
    }

    bool HasNaNs() {
        return isNaN(x) || isNaN(y) || isNaN(z);
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

template <typename T> inline Vector3<T> Cross(const Vector3<T> &v, const Normal3<T> &n) {
    double vx = v.x, vy = v.y, vz = v.z;
    double nx = n.x, ny = n.y, nz = n.z;
    return Vector3<T>(
        (vy * nz) - (vz * ny),
        (vz * nx) - (vx * nz),
        (vx * ny) - (vy * nx)
    );
}

template <typename T> inline Vector3<T> Cross(const Normal3<T> &v, const Vector3<T> &n) {
    double vx = v.x, vy = v.y, vz = v.z;
    double nx = n.x, ny = n.y, nz = n.z;
    return Vector3<T>(
        (vy * nz) - (vz * ny),
        (vz * nx) - (vx * nz),
        (vx * ny) - (vy * nx)
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

    Point2<T> &operator=(const Point2<T> &p) {
        x = p.x;
        y = p.y;
        return *this;
    }

    bool operator==(const Point2<T> &p) const {
        return x == p.x && y == p.y;
    }

    bool operator!=(const Point2<T> &p) const {
        return x != p.x || y != p.y;
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
        return Point2<T>(x+p.x, y+p.y);
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

    Point2<T> operator/(T f) const {
        Assert(f != 0);
        // Multiplication is faster than division. Compilers are generally restricted
        // from performing transformations of this type. So they wouldn't happen automatically.
        Float reciprocal = (Float)1/f;
        return Point2<T>(x*reciprocal, y*reciprocal);
    }

    T operator[](int i) const {
        Assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    T &operator[](int i) {
        Assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    bool HasNaNs() const {
        return isNaN(x) || isNaN(y);
    }
};

typedef Point2<Float> Point2f;
typedef Point2<int> Point2i;

template <typename T> inline Float Distance(const Point2<T> &p1, const Point2<T> &p2) {
    // Guaranteed to be positive.
    return (p1-p2).Length();
}

template <typename T> inline Float DistanceSquared(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1-p2).LengthSquared();
}

template <typename T> inline Point2<T> operator*(T s, const Point2<T> &p) {
    return p * s;
}

template <typename T, typename U> inline Point2<T> operator*(U f, const Point2<T> &p) {
    Assert(!p.HasNaNs());
    return p * f;
}

// Linear interpolation (0 <= t <= 1) and extrapolation (t < 0, t > 1).
template <typename T> inline Point2<T> Lerp(Float t, const Point2<T> &p0, const Point2<T> &p1) {
    return (1-t)*p0 + t*p1;
}

// Component-wise minimum.
template <typename T> inline Point2<T> Min(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(
        std::min(p1.x, p2.x),
        std::min(p1.y, p2.y)
    );
}

// Component-wise maximum.
template <typename T> inline Point2<T> Max(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(
        std::max(p1.x, p2.x),
        std::max(p1.y, p2.y)
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

    Point3<T> &operator=(const Point3<T> &p) {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    bool operator==(const Point3<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const Point3<T> &p) const {
        return x != p.x || y != p.y || y != p.z;
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

    template <typename U> Point3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }

    Point3<T> operator/(T f) const {
        Assert(f != 0);
        // Multiplication is faster than division. Compilers are generally restricted
        // from performing transformations of this type. So they wouldn't happen automatically.
        Float reciprocal = (Float)1/f;
        return Point3<T>(x*reciprocal, y*reciprocal, z*reciprocal);
    }

    T operator[](int i) const {
        Assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    T &operator[](int i) {
        Assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    bool HasNaNs() {
        return isNaN(x) || isNaN(y);
    }
};

typedef Point3<Float> Point3f;
typedef Point3<int> Point3i;

template <typename T, typename U> inline Point3<T> operator*(U f, const Point3<T> &p) {
    return p * f;
}

template <typename T> inline Float Distance(const Point3<T> &p1, const Point3<T> &p2) {
    // Guaranteed to be positive.
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
        Assert(!HasNaNs());
    }

    explicit Normal3<T>(const Vector3<T> &v)
        : x(v.x), y(v.y), z(v.z) 
    {
        Assert(!HasNaNs());
    }

    bool operator==(const Normal3<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const Normal3<T> &p) const {
        return x != p.x || y != p.y || y != p.z;
    }

    template <typename U> Normal3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
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

    template <typename U> Normal3<T> operator/(U f) const {
        Float inv = (Float) 1 / f;
        return Normal3<T>(x * inv, y * inv, z * inv);
    }

    Float LengthSquared() const {
        return x*x + y*y + z*z;
    }

    Float Length() const {
        return std::sqrt(LengthSquared());
    }

    bool HasNaNs() {
        return isNaN(x) || isNaN(y) || isNaN(z);
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

template <typename T> inline T Dot(const Normal3<T> &n1, const Normal3<T> &n2) {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T> inline T AbsDot(const Normal3<T> &n1, const Normal3<T> &n2) {
    return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

template <typename T> inline T Dot(const Normal3<T> &n, const Vector3<T> &v) {
    return n.x*v.x + n.y*v.y + n.z*v.z;
}

template <typename T> inline T AbsDot(const Normal3<T> &n, const Vector3<T> &v) {
    return std::abs(Dot(n, v));
}

template <typename T> inline T Dot(const Vector3<T> &v, const Normal3<T> &n) {
    return v.x * n.x + v.y * n.y + v.z * n.z;
}

template <typename T> inline T AbsDot(const Vector3<T> &v, const Normal3<T> &n) {
    return std::abs(Dot(v, n));
}

// Flip a normal so that it lies in the same hemisphere as the vector.
template <typename T> inline Normal3<T> Faceforward(const Normal3<T> &n, const Vector3<T> &v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T> inline Normal3<T> Faceforward(const Normal3<T> &n, const Normal3<T> &n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T> inline Vector3<T> Faceforward(const Vector3<T> &v, const Vector3<T> &v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T> inline Vector3<T> Faceforward(const Vector3<T> &v, const Normal3<T> &n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template <typename T> Normal3<T> Abs(const Normal3<T> &v) {
    return Normal3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

class Ray {
public:
    // Origin.
    Point3f o;

    // Direction.
    Vector3f d;

    // Restricts the ray to a segment of points [0, r(tMax)].
    mutable Float tMax;

    Float time;

    // The medium containing the ray's origin.
    const Medium *medium;

    // Point3F and Vector3f default constructors set the origin and direction to 0.
    Ray() : tMax(Infinity), time(0.f), medium(nullptr) {}

    Ray(
        const Point3f &o,
        const Vector3f &d,
        Float tMax = Infinity,
        Float time = 0.f,
        const Medium *medium = nullptr
    ) : o(o), d(d), tMax(tMax), time(time), medium(medium) {}

    Ray(const Ray& r) : o(r.o), d(r.d), tMax(r.tMax), time(r.time), medium(r.medium) {}

    // Overloaded function call operator returns the point along the ray at t.
    // Ex. Ray r(...); Float t = 1.3; Point3f p = r(t); 
    Point3f operator()(Float t) const {
        return o + d*t;
    }
};

class RayDifferential : public Ray {
public:
    bool hasDifferentials;

    // Auxiliary ray offset in the x direction.
    Point3f rxOrigin;
    Vector3f rxDirection;

    // Auxiliary ray offset in the y direction.
    Point3f ryOrigin;
    Vector3f ryDirection;

    RayDifferential() {
        hasDifferentials = false;
    }

    RayDifferential(
        const Point3f &o,
        const Vector3f &d,
        Float tMax = Infinity,
        Float time = 0.f,
        const Medium *medium = nullptr
    ) : Ray(o, d, tMax, time, medium) {
        hasDifferentials = false;
    }

    RayDifferential(const Ray &ray) : Ray(ray) {
        hasDifferentials = false;
    }

    // Scales the distance between the main ray and the x and y ray differentials.
    // A RayDifferential is used to sample the 2 neighboring pixels of the pixel sampled
    // by the main ray. By default, the differentials and the main ray are 1 pixel apart,
    // but this distance may be changed by this method: in the case when multiple samples
    // of a pixel are taken (supersampling), we want the differentials to be closer to the
    // main ray; the more samples per pixel are taken, the shorter the distance should be.
    void ScaleDifferentials(Float s) {
        rxOrigin = o + (rxOrigin-o)*s;
        ryOrigin = o + (ryOrigin-o)*s;
        rxDirection = d + (rxDirection-d)*s;
        ryDirection = d + (ryDirection-d)*s;
    }
};

template <typename T> class Bounds2 {
public:
    // Invariant: pMin.x <= pMax.x AND pMin.y <= pMax.y.
    Point2<T> pMin;
    Point2<T> pMax;

    Bounds2() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        // Uninitialized bounding box violates invariant: 
        // pMin.x > pMax.x, pMin.y > pMax.y.
        // The union or intersection of this BB (with a valid BB) will see it empty, which would
        // yield the correct results.
        pMin = Point2<T>(maxNum, maxNum);
        pMax = Point2<T>(minNum, minNum);
    }

    // Enclose a single point.
    Bounds2(const Point2<T> &p) : pMin(p), pMax(p) {}

    // p1 and p2 may indeed be the opposite corner points of a box, but they may not have been
    // passed in the order that satisfies the pMin, pMax invariant.
    Bounds2(const Point2<T> &p1, const Point2<T> &p2)
        : pMin(Point2<T>(Min(p1, p2))),
          pMax(Point2<T>(Max(p1, p2)))
    {}

    // Cast between Bounds2 types.
    template <typename U> explicit operator Bounds2<U>() const {
        return Bounds2<U>((Point2<U>) pMin, (Point2<U>) pMax);
    }

    Point2<T> operator[](int i) const {
        Assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }

    Point2<T> &operator[](int i) {
        Assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }

    // Coordinates of the input corner, where pMin is corner 0 and pMax is corner 3.
    Point2<T> Corner(int corner) const {
        // If corner = 3 = 0b11, returned point is (pMax.x, pMax.y).
        // If corner = 1 = 0b01, returned point is (pMin.x, pMax.y).
        return Point2<T>(
            (*this)[(corner & 1)].x,
            (*this)[(corner & 2) ? 1: 0].y
        );
    }

    Vector2<T> Diagonal() const {
        return pMax - pMin;
    }

    T Area() const {
        Vector2<T> d = Diagonal();
        return d.x*d.y;
    }

    // Returns the index of the longest axis.
    int MaximumExtent() const {
        Vector2<T> = Diagonal();
        if (d.x > d.y) {
            return 0;
        }
        return 1;
    }

    // The interpolation parameter t allows for a different parameter for each dimension. 
    Point2<T> Lerp(const Point2f &t) const {
        return Point2<T>(
            ::Lerp(t.x, pMin.x, pMax.x),
            ::Lerp(t.y, pMin.y, pMax.y)
        );
    }

    // When p lies inside the box, its offset is its position relative to the pMin and pMax
    // corners of the box along each of the axes. The offset is of the form 
    // ([0.0, 1.0], [0.0, 1.0]), the offset of a point exactly at pMin being (0.0, 0.0) and
    // the offset of one exactly at pMax being (1.0, 1.0). For example, a point with offset
    // (0.2, 0.7) is closer to pMin along the x-axis and closer to pMax along the y-axis.
    Vector2<T> Offset(const Point2<T> &p) const {
        Vector2<T> o = p - pMin;
        if (pMax.x > pMin.x) {
            o.x /= pMax.x - pMin.x;
        }
        if (pMax.y > pMin.y) {
            o.y /= pMax.y - pMin.y;
        }
        return o;
    }

    void BoundingCircle(Point2<T> *center, Float *radius) const {
        *center = (pMax + pMin) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }
};

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;

// Construct a new extended box that contains point p.
template <typename T> inline Bounds2<T> Union(const Bounds2<T> &b, const Point2<T> &p) {
    return Bounds2<T>(
        Point2<T>(std::min(b.pMin.x, p.x), std::min(b.pMin.y, p.y)),
        Point2<T>(std::max(b.pMax.x, p.x), std::max(b.pMax.y, p.y))
    );
}

// Construct a new extended box that bounds the space of the 2 boxes.
template <typename T> inline Bounds2<T> Union(const Bounds2<T> &b1, const Bounds2<T> &b2) {
    return Bounds2<T>(
        Point2<T>(std::min(b1.pMin.x, b2.pMin.x), std::min(b1.pMin.y, b2.pMin.y)),
        Point2<T>(std::max(b1.pMax.x, b2.pMax.x), std::max(b1.pMax.y, b2.pMax.y))
    );
}

// Construct a new box that bounds the space of the intersection of the 2 boxes.
template <typename T> inline Bounds2<T> Intersect(const Bounds2<T> &b1, const Bounds2<T> &b2) {
    return Bounds2<T>(
        Point2<T>(std::max(b1.pMin.x, b2.pMin.x), std::max(b1.pMin.y, b2.pMin.y)),
        Point2<T>(std::min(b1.pMax.x, b2.pMax.x), std::min(b1.pMax.y, b2.pMax.y))
    );
}

template <typename T> inline bool Overlaps(const Bounds2<T> &b1, const Bounds2<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    return x && y;
}

template <typename T> inline bool Inside(const Point2<T> &p, const Bounds2<T> &b) {
    return p.x >= b.pMin.x && p.x <= b.pMax.x
        && p.y >= b.pMin.y && p.y <= b.pMax.y;
}

// Doesn't consider points on the upper boundary to be inside the bounds. Mostly useful
// for interger-typed bounds.
template <typename T> inline bool InsideExclusive(const Point2<T> &p, const Bounds2<T> &b) {
    return p.x >= b.pMin.x && p.x < b.pMax.x
        && p.y >= b.pMin.y && p.y < b.pMax.y;
}

template <typename T, typename U> inline Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
    return Bounds2<T>(
        b.pMin - Vector2<T>(delta, delta),
        b.pMax + Vector2<T>(delta, delta)
    );
}

class Bounds2iIterator : public std::forward_iterator_tag {
private:
    Point2i p;

    const Bounds2i *bounds;

    void advance() {
        ++p.x;
        if (p.x == bounds->pMax.x) {
            p.x = bounds->pMin.x;
            ++p.y;
        }
    }

public:
    Bounds2iIterator(const Bounds2i &b, const Point2i &p) : bounds(&b), p(p) {}

    // Postfix increment.
    Bounds2iIterator operator++() {
        advance();
        return *this;
    }

    // Prefix increment.
    Bounds2iIterator operator++(int) {
        Bounds2iIterator old = *this;
        advance();
        return old;
    }

    bool operator==(const Bounds2iIterator &bi) const {
        return p == bi.p && bounds == bi.bounds;
    }

    bool operator!=(const Bounds2iIterator &bi) const {
        return p != bi.p || bounds != bi.bounds;
    }

    Point2i operator*() const {
        return p;
    }
};

inline Bounds2iIterator begin(const Bounds2i &b) {
    return Bounds2iIterator(b, b.pMin);
}

inline Bounds2iIterator end(const Bounds2i &b) {
    Point2i pEnd(b.pMin.x, b.pMax.y);

    if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y) {
        // Degenerate.
        pEnd = b.pMin;
    }

    return Bounds2iIterator(b, pEnd);
}

template <typename T> class Bounds3 {
public:
    // Invariant: pMin.x <= pMax.x AND pMin.y <= pMax.y AND pMin.z <= pMax.z.
    Point3<T> pMin;
    Point3<T> pMax;

    Bounds3() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        // Uninitialized bounding box violates invariant: 
        // pMin.x > pMax.x, pMin.y > pMax.y, pMin.z <= pMax.z.
        // The union or intersection of this BB (with a valid BB) will see it empty, which would
        // yield the correct results.
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }

    // Enclose a single point.
    Bounds3(const Point3<T> &p) : pMin(p), pMax(p) {}

    // p1 and p2 may indeed be the opposite corner points of a box, but they may not have been
    // passed in the order that satisfies the pMin, pMax invariant.
    Bounds3(const Point3<T> &p1, const Point3<T> &p2)
        : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
          pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z))
    {}

    Point3<T> operator[](int i) const {
        Assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }

    Point3<T> &operator[](int i) {
        Assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }

    // Coordinates of the input corner, where pMin is corner 0 and pMax is corner 7.
    Point3<T> Corner(int corner) const {
        // If corner = 6 = 0b110, returned point is (pMin.x, pMax.y, pMax.z).
        // If corner = 4 = 0b100, returned point is (pMin.x, pMin.y, pMax.z).
        return Point3<T>(
            (*this)[(corner & 1)].x,
            (*this)[(corner & 2) ? 1: 0].y,
            (*this)[(corner & 4) ? 1 : 0].z
        );
    }

    Vector3<T> Diagonal() const {
        return pMax - pMin;
    }

    T SurfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x*d.y + d.x*d.z + d.y*d.z);
    }

    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }

    // Returns the index of the longest axis. Useful when deciding which axis to
    // subdivide when building some of the ray-intersection acceleration structures.
    int MaximumExtent() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z) {
            return 0;
        } else if (d.y > d.z) {
            return 1;
        }
        return 2;
    }

    // The interpolation parameter t allows for a different parameter for each dimension. 
    Point3<T> Lerp(const Point3f &t) const {
        return Point3<T>(
            ::Lerp(t.x, pMin.x, pMax.x),
            ::Lerp(t.y, pMin.y, pMax.y),
            ::Lerp(t.z, pMin.z, pMax.z)
        );
    }

    // When p lies inside the box, its offset is its position relative to the pMin and pMax
    // corners of the box along each of the axes. The offset is of the form 
    // ([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]), the offset of a point exactly at pMin being 
    // (0.0, 0.0, 0.0) and the offset of one exactly at pMax being (1.0, 1.0, 1.0). For 
    // example, a point with offset (0.2, 0.7, 0.5) is closer to pMin along the x-axis,
    // closer to pMax along the y-axis, and halfway between pMin and pMax along the z-axis.
    Vector3<T> Offset(const Point3<T> &p) const {
        Vector3<T> o = p - pMin;
        if (pMax.x > pMin.x) {
            o.x /= pMax.x - pMin.x;
        }
        if (pMax.y > pMin.y) {
            o.y /= pMax.y - pMin.y;
        }
        if (pMax.z > pMin.z) {
            o.z /= pMax.z - pMin.z;
        }
        return o;
    }

    void BoundingSphere(Point3<T> *center, Float *radius) const {
        *center = (pMax + pMin) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }

    // P is for "predicate".
    bool IntersectP(const Ray &ray, Float *hitt0 = nullptr, Float *hitt1 = nullptr) const;
    bool IntersectP(const Ray &ray, const Vector3f &reciprocalDir, const int dirIsNeg[3]) const;
};

typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

template <typename T> inline bool Bounds3<T>::IntersectP(
    const Ray &ray,
    Float *hitt0,
    Float *hitt1
) const {
    // The parametric interval [t0, t1] shrinks with every slab intersection that is found.
    Float t0 = 0;
    Float t1 = ray.tMax;

    // Iteration 1 tests the intersection of the ray with the slab between the planes that
    // are perpendicular to the X axis. Iteration 2, Y axis. Iteration 3, Z axis.
    for (int i=0; i<3; ++i) {
        // Update interval for ith bounding box slab.

        // Avoid division later on. If component is 0, division by 0 of any nonzero floating-
        // point number results in a special +infinity or -infinity value in IEEE floating-
        // point arithmetic. The algorithm is supposed to be correct in the face of such case.
        Float reciprocalDirComponent = 1 / ray.d[i];
        
        // Compute values of the ray's parameter for which the ray satisfies the implicit
        // equation of the plane of. The ray intersects the slab at points (r.o + tNear*r.d)
        // and (r.o + tFar*r.d), if at all.
        Float tNear = (pMin[i] - ray.o[i]) *  reciprocalDirComponent;
        Float tFar = (pMax[i] - ray.o[i]) * reciprocalDirComponent;
        if (tNear > tFar) std::swap(tNear, tFar);

        // Update tFar to ensure robust ray-bounds intersection (see 3.9.2).
        tFar *= 1 + 2 * gamma(3);

        // Update (shrink) the parametric interval from slab intersection t values: 
        // [t0, t1] = [tNear, tFar].
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1) {
            // The parametric interval is empty; it shrunk until the endpoints crossed each other.
            return false;
        }
    }

    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}

// Optimized overload. 15% performance improvement approximately.
template <typename T> inline bool Bounds3<T>::IntersectP(
    const Ray &ray, 
    const Vector3f &reciprocalDir, 
    const int dirIsNeg[3]
) const {
    // TODO: implement.
    return false;
}

// Construct a new extended box that contains point p.
template <typename T> inline Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
    return Bounds3<T>(
        Point3<T>(std::min(b.pMin.x, p.x), std::min(b.pMin.y, p.y), std::min(b.pMin.z, p.z)),
        Point3<T>(std::max(b.pMax.x, p.x), std::max(b.pMax.y, p.y), std::max(b.pMax.z, p.z))
    );
}

// Construct a new extended box that bounds the space of the 2 boxes.
template <typename T> inline Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    return Bounds3<T>(
        Point3<T>(std::min(b1.pMin.x, b2.pMin.x), std::min(b1.pMin.y, b2.pMin.y), std::min(b1.pMin.z, b2.pMin.z)),
        Point3<T>(std::max(b1.pMax.x, b2.pMax.x), std::max(b1.pMax.y, b2.pMax.y), std::max(b1.pMax.z, b2.pMax.z))
    );
}

// Construct a new box that bounds the space of the intersection of the 2 boxes.
template <typename T> inline Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    return Bounds3<T>(
        Point3<T>(std::max(b1.pMin.x, b2.pMin.x), std::max(b1.pMin.y, b2.pMin.y), std::max(b1.pMin.z, b2.pMin.z)),
        Point3<T>(std::min(b1.pMax.x, b2.pMax.x), std::min(b1.pMax.y, b2.pMax.y), std::min(b1.pMax.z, b2.pMax.z))
    );
}

template <typename T> inline bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return x && y && z;
}

template <typename T> inline bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
    return p.x >= b.pMin.x && p.x <= b.pMax.x
        && p.y >= b.pMin.y && p.y <= b.pMax.y
        && p.z >= b.pMin.z && p.z <= b.pMax.z;
}

// Doesn't consider points on the upper boundary to be inside the bounds. Mostly useful
// for interger-typed bounds.
template <typename T> inline bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
    return p.x >= b.pMin.x && p.x < b.pMax.x
        && p.y >= b.pMin.y && p.y < b.pMax.y
        && p.z >= b.pMin.z && p.z < b.pMax.z;
}

template <typename T, typename U> inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    return Bounds3<T>(
        b.pMin - Vector3<T>(delta, delta, delta),
        b.pMax + Vector3<T>(delta, delta, delta)
    );
}

// Converts a spherical coordinate to a rectangular coordinate in a standard coordinate system.
// Theta is measured from the z axis. Phi is measured about the z axis from the x axis.
inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

// Converts a spherical coordinate to a rectangular coordinate in the coordinate frame given
// by the x, y, and z basis vectors. Theta is measured from the z axis. Phi is measured about
// the z axis from the x axis.
inline Vector3f SphericalDirection(
    Float sinTheta,
    Float cosTheta,
    Float phi,
    const Vector3f &x,
    const Vector3f &y,
    const Vector3f &z
) {
    return (sinTheta * std::cos(phi) * x) + (sinTheta * std::sin(phi) * y) + (cosTheta * z);
}

// Converts a rectangular coordinate to a spherical theta angle (measured from the z axis).
// z = cos(theta), so theta = arccos(z).
inline Float SphericalTheta(const Vector3f &v) {
    // v is assumed to be normalized.
    return std::acos(Clamp(v.z, -1, 1));
}

// Converts a rectangular coordinate to a spherical phi angle (measured about the z axis from
// the x axis).
inline Float SphericalPhi(const Vector3f &v) {
    Float phi = std::atan2(v.y, v.x);
    return (phi < 0) ? (phi + 2*Pi) : phi;
}

inline Point3f OffsetRayOrigin(
    const Point3f &p,
    const Vector3f &pError,
    const Normal3f &n,
    const Vector3f &w
) {
    Float d = Dot(Abs(n), pError);
    Vector3f offset = d * Vector3f(n);
    if (Dot(w, n) < 0) offset = -offset;
    Point3f po = p + offset;
    // Round offset point _po_ away from _p_.
    for (int i = 0; i < 3; ++i) {
        if (offset[i] > 0) {
            po[i] = NextFloatUp(po[i]);
        }
        else if (offset[i] < 0) {
            po[i] = NextFloatDown(po[i]);
        }
    }
    return po;
}

#endif // CPBRT_CORE_GEOMETRY_H