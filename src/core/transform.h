#ifndef CPBRT_CORE_TRANSFORM_H
#define CPBRT_CORE_TRANSFORM_H

#include "cpbrt.h"
#include "geometry.h"

class Matrix4x4 {
public:
    // Row-major order.
    Float m[4][4];

    Matrix4x4() {
        // Identity matrix.
        m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
        m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
        m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
        m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
    }

    Matrix4x4(const Float m[4][4]) {
        memcpy(this->m, m, 16*sizeof(Float));
    }

    Matrix4x4(
        Float t00, Float t01, Float t02, Float t03,
        Float t10, Float t11, Float t12, Float t13,
        Float t20, Float t21, Float t22, Float t23,
        Float t30, Float t31, Float t32, Float t33
    ) {
        m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
        m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
        m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
        m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;
    }

    bool operator==(const Matrix4x4 &mat) const {
        return 
            m[0][0] == mat.m[0][0] && m[0][1] == mat.m[0][1] && m[0][2] == mat.m[0][2] && m[0][3] == mat.m[0][3] &&
            m[1][0] == mat.m[1][0] && m[1][1] == mat.m[1][1] && m[1][2] == mat.m[1][2] && m[1][3] == mat.m[1][3] &&
            m[2][0] == mat.m[2][0] && m[2][1] == mat.m[2][1] && m[2][2] == mat.m[2][2] && m[2][3] == mat.m[2][3] &&
            m[3][0] == mat.m[3][0] && m[3][1] == mat.m[3][1] && m[3][2] == mat.m[3][2] && m[3][3] == mat.m[3][3];
    }

    bool operator!=(const Matrix4x4 &mat) const {
        return 
            m[0][0] != mat.m[0][0] || m[0][1] != mat.m[0][1] || m[0][2] != mat.m[0][2] || m[0][3] != mat.m[0][3] ||
            m[1][0] != mat.m[1][0] || m[1][1] != mat.m[1][1] || m[1][2] != mat.m[1][2] || m[1][3] != mat.m[1][3] ||
            m[2][0] != mat.m[2][0] || m[2][1] != mat.m[2][1] || m[2][2] != mat.m[2][2] || m[2][3] != mat.m[2][3] ||
            m[3][0] != mat.m[3][0] || m[3][1] != mat.m[3][1] || m[3][2] != mat.m[3][2] || m[3][3] != mat.m[3][3];
    }

    Matrix4x4 operator*(const Matrix4x4 &mat) const;
    friend Matrix4x4 Mul(const Matrix4x4 &mat1, const Matrix4x4 &mat2);

    Matrix4x4 Transpose() const;
    friend Matrix4x4 Transpose(const Matrix4x4 &mat);

    Matrix4x4 Inverse() const;
    friend Matrix4x4 Inverse(const Matrix4x4 &mat);
};

class Transform {
private:
    // Row-major order.
    Matrix4x4 m;

    // Inverse of m.
    Matrix4x4 mInv;

public:
    // Identity transformation (matrix defaults to identity).
    // The identity matrix is its own inverse.
    Transform() {}

    Transform(const Float m[4][4]) {
        this->m = Matrix4x4(m);

        mInv = this->m.Inverse();
    }

    Transform(const Matrix4x4 &mat) : m(mat), mInv(mat.Inverse()) {}

    // Let the caller supply the inverse. If the caller passes an inverse written
    // by hand, the saved expense of computing it is significant.
    Transform(const Matrix4x4 &mat, const Matrix4x4 &matInv) : m(mat), mInv(matInv) {}

    Transform Inverse() const {
        return Transform(mInv, m);
    }

    friend Transform Inverse(const Transform &t) {
        return t.Inverse();
    }

    Transform Transpose() const {
        return Transform(m.Transpose(), mInv.Transpose());
    }

    const Matrix4x4 &GetMatrix() const { 
        return m;
    }

    bool IsIdentity() const {
        return m == Matrix4x4();
    }

    bool operator==(const Transform &t) const {
        // The inverse of a matrix is unique; the second condition is unnecessary
        // unless the inverse wasn't computed with Matrix4x4::Inverse() and happens
        // to be wrong.
        return m == t.m && mInv == t.mInv;
    }

    bool operator!=(const Transform &t) const {
        return m != t.m || mInv != t.mInv;
    }

    template <typename T> Point3<T> operator()(const Point3<T> &p) const;

    template <typename T> inline Point3<T> operator()(const Point3<T> &pt, Vector3<T> *absError) const;

    template <typename T> inline Point3<T> operator()(const Point3<T> &p, const Vector3<T> &pError, Vector3<T> *pTransError) const;

    template <typename T> Vector3<T> operator()(const Vector3<T> &v) const;

    template <typename T> inline Vector3<T> operator()(const Vector3<T> &v, Vector3<T> *vTransError) const;

    template <typename T> inline Vector3<T> operator()(const Vector3<T> &v, const Vector3<T> &vError, Vector3<T> *vTransError) const;

    template <typename T> Normal3<T> operator()(const Normal3<T> &n) const;

    inline Ray operator()(const Ray &r) const;

    inline Ray operator()(const Ray &r, Vector3f *oError, Vector3f *dError) const;

    inline Ray operator()(
        const Ray &r,
        const Vector3f &oErrorIn, const Vector3f &dErrorIn,
        Vector3f *oErrorOut, Vector3f *dErrorOut
    ) const;

    inline RayDifferential operator()(const RayDifferential &dr) const;

    template <typename T> Bounds3<T> operator()(const Bounds3<T> &aabb) const;

    SurfaceInteraction operator()(const SurfaceInteraction &si) const;

    Transform operator*(const Transform &t) const;

    Transform Translate(const Vector3f &delta) const;

    Transform Scale(Float sx, Float sy, Float sz) const;

    // Determines if the transformation has a scaling factor on any dimension,
    // (only when its matrix is the identity does it not have a scaling factor).
    bool HasScale() const {
        Float la2 = (*this)(Vector3f(1, 0, 0)).LengthSquared();
        Float lb2 = (*this)(Vector3f(0, 1, 0)).LengthSquared();
        Float lc2 = (*this)(Vector3f(0, 0, 1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef  NOT_ONE
    }

    Transform RotateX(Float theta) const;

    Transform RotateY(Float theta) const;

    Transform RotateZ(Float theta) const;

    Transform Rotate(Float theta, const Vector3f & axis) const;

    // Transforms a point in camera space to world space. The camera space's frame is
    // (pos, ?, up, look), where the iHat basis vector is implicitly defined by the others. 
    Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up) const;

    bool SwapsHandedness() const;
};

template <typename T> inline Point3<T> Transform::operator()(const Point3<T> &p) const {
    T homogeneousX = (m.m[0][0] * p.x) + (m.m[0][1] * p.y) + (m.m[0][2] * p.z) + m.m[0][3];
    T homogeneousY = (m.m[1][0] * p.x) + (m.m[1][1] * p.y) + (m.m[1][2] * p.z) + m.m[1][3];
    T homogeneousZ = (m.m[2][0] * p.x) + (m.m[2][1] * p.y) + (m.m[2][2] * p.z) + m.m[2][3];
    T weight       = (m.m[3][0] * p.x) + (m.m[3][1] * p.y) + (m.m[3][2] * p.z) + m.m[3][3];

   if (weight != 1.0f) {
       return Point3<T>(
           homogeneousX, 
           homogeneousY,
           homogeneousZ
       ) / weight;
   }
   return Point3<T>(homogeneousX, homogeneousY, homogeneousZ);
};

template <typename T> inline Point3<T> Transform::operator()(const Point3<T> &p, Vector3<T> *pError) const {
    T x = p.x, y = p.y, z = p.z;
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    T xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                 std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    T yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                 std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    T zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                 std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    *pError = gamma(3) * Vector3<T>(xAbsSum, yAbsSum, zAbsSum);

    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template <typename T> inline Point3<T> Transform::operator()(
    const Point3<T> &pt,
    const Vector3<T> &ptError,
    Vector3<T> *absError
) const {
    T x = pt.x, y = pt.y, z = pt.z;
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * ptError.x + std::abs(m.m[0][1]) * ptError.y +
             std::abs(m.m[0][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                    std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * ptError.x + std::abs(m.m[1][1]) * ptError.y +
             std::abs(m.m[1][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                    std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * ptError.x + std::abs(m.m[2][1]) * ptError.y +
             std::abs(m.m[2][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                    std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    if (wp == 1.)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template <typename T> inline Vector3<T> Transform::operator()(const Vector3<T> &v) const {
    // No need to use the homogeneous representation of vectors to transform them,
    // because their weight, 0, causes entry (mv)30 to be 0.
    return Vector3<T>(
        (m.m[0][0] * v.x) + (m.m[0][1] * v.y) + (m.m[0][2] * v.z),
        (m.m[1][0] * v.x) + (m.m[1][1] * v.y) + (m.m[1][2] * v.z),
        (m.m[2][0] * v.x) + (m.m[2][1] * v.y) + (m.m[2][2] * v.z)
    );
};

template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v, Vector3<T> *absError) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x = gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) + std::abs(m.m[0][2] * v.z));
    absError->y = gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) + std::abs(m.m[1][2] * v.z));
    absError->z = gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) + std::abs(m.m[2][2] * v.z));
    return Vector3<T>(
        m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z
    );
}

template <typename T> inline Vector3<T> Transform::operator()(
    const Vector3<T> &v,
    const Vector3<T> &vError,
    Vector3<T> *absError
) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * vError.x + std::abs(m.m[0][1]) * vError.y +
             std::abs(m.m[0][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
                    std::abs(m.m[0][2] * v.z));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * vError.x + std::abs(m.m[1][1]) * vError.y +
             std::abs(m.m[1][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
                    std::abs(m.m[1][2] * v.z));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * vError.x + std::abs(m.m[2][1]) * vError.y +
             std::abs(m.m[2][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
                    std::abs(m.m[2][2] * v.z));
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T> inline Normal3<T> Transform::operator()(const Normal3<T> &n) const {
    // The transform doesn't maintain the orthogonality of the normal and the surface 
    // (the tangent to the surface). The inverse transpose does.
    T x = n.x, y = n.y, z = n.z;
    return Normal3<T>(
        mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
        mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
        mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z
    );
}

inline Ray Transform::operator()(const Ray &r) const {
    Vector3f oError;
    Point3f o = (*this)(r.o, &oError);
    Vector3f d = (*this)(r.d);
    Float lengthSquared = d.LengthSquared();
    Float tMax = r.tMax;
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), oError) / lengthSquared;
        o += d * dt;
        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

inline RayDifferential Transform::operator()(const RayDifferential &rd) const {
    Ray transformedRay = (*this)((Ray) rd);
    RayDifferential transformedDiff = RayDifferential(transformedRay);
    transformedDiff.hasDifferentials = rd.hasDifferentials;
    transformedDiff.rxOrigin = (*this)(rd.rxOrigin);
    transformedDiff.rxDirection = (*this)(rd.rxDirection);
    transformedDiff.ryOrigin = (*this)(rd.ryOrigin);
    transformedDiff.ryDirection = (*this)(rd.ryDirection);
    return transformedDiff;
}

inline Ray Transform::operator()(const Ray &r, Vector3f *oError, Vector3f *dError) const {
    Point3f o = (*this)(r.o, oError);
    Vector3f d = (*this)(r.d, dError);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oError) / lengthSquared;
        o += d * dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

inline Ray Transform::operator()(
    const Ray &r,
    const Vector3f &oErrorIn,
    const Vector3f &dErrorIn,
    Vector3f *oErrorOut,
    Vector3f *dErrorOut
) const {

    Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
    Vector3f d = (*this)(r.d, dErrorIn, dErrorOut);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oErrorOut) / lengthSquared;
        o += d * dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

template <typename T> inline Bounds3<T> Transform::operator()(const Bounds3<T> &aabb) const {
    Bounds3<T> transformedAABB((*this)(aabb.Corner(0)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(1)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(2)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(3)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(4)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(5)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(6)));
    transformedAABB = Union(transformedAABB, (*this)(aabb.Corner(7)));
    return transformedAABB;
}

Transform Translate(const Vector3f &delta);

Transform Scale(Float x, Float y, Float z);

Transform RotateX(Float theta);

Transform RotateY(Float theta);

Transform RotateZ(Float theta);

Transform Rotate(Float theta, const Vector3f &axis);

Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up);

// Solves the Ax=b linear system using Cramer's rule.
bool SolveLinearSystem2x2(const Float A[2][2], const Float b[2], Float *x0, Float *x1);

#endif // CPBRT_CORE_TRANSFORM_H