#ifndef CPBRT_CORE_TRANSFORM_H
#define CPBRT_CORE_TRANSFORM_H

#include "cpbrt.h"
#include "geometry.h"
#include "interaction.h"

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

    template <typename T> Vector3<T> operator()(const Vector3<T> &v) const;

    template <typename T> Normal3<T> operator()(const Normal3<T> &n) const;

    Ray operator()(const Ray &r) const;

    RayDifferential operator()(const RayDifferential &dr) const;

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

// Solves the Ax=b linear system using Cramer's rule.
bool SolveLinearSystem2x2(const Float A[2][2], const Float b[2], Float *x0, Float *x1) {
    Float detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (std::abs(detA) < 1e-10f) {
        return false;
    }

    // A0b is the 2x2 matrix formed by replacing column 0 of A with b.
    Float detA0b = b[0]*A[1][1] - A[0][1]*b[1];
    *x0 = detA0b / detA;

    // A1b is the 2x2 matrix formed by replacing column 1 of A with b.
    Float detA1b = A[0][0]*b[1] - b[0]*A[1][0];
    *x1 = detA1b / detA;

    if (std::isnan(*x1) || std::isnan(*x1)) {
        return false;
    }
    return true;
}

#endif // CPBRT_CORE_TRANSFORM_H