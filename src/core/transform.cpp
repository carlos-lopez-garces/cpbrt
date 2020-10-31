#include "transform.h"

// Matrix4x4 method definitions.

Matrix4x4 Matrix4x4::Transpose() const {
    return Matrix4x4(
        m[0][0], m[1][0], m[2][0], m[3][0],
        m[0][1], m[1][1], m[2][1], m[3][1],
        m[0][2], m[1][2], m[2][2], m[3][2],
        m[0][3], m[1][3], m[2][3], m[3][3]
    );
}

Matrix4x4 Transpose(const Matrix4x4 &mat) {
    return Matrix4x4(
        mat.m[0][0], mat.m[1][0], mat.m[2][0], mat.m[3][0],
        mat.m[0][1], mat.m[1][1], mat.m[2][1], mat.m[3][1],
        mat.m[0][2], mat.m[1][2], mat.m[2][2], mat.m[3][2],
        mat.m[0][3], mat.m[1][3], mat.m[2][3], mat.m[3][3]
    );
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4 &mat) const {
    Matrix4x4 product;

    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            product.m[i][j] = 
                m[i][0] * mat.m[0][j] 
                + m[i][1] * mat.m[1][j] 
                + m[i][2] * mat.m[2][j]
                + m[i][3] * mat.m[3][j];
        }
    }

    return product;
}

Matrix4x4 Mul(const Matrix4x4 &mat1, const Matrix4x4 &mat2) {
    return mat1 * mat2;
}

// Transform method definitions.

template <typename T> inline Point3<T> Transform::operator()(const Point3<T> &p) const {
    T homogeneousX = (m.m[0][0] * p.x) + (m.m[0][1] * p.y) + (m.m[0][2] * p.z) + m.m[0][3];
    T homogeneousY = (m.m[1][0] * p.x) + (m.m[1][1] * p.y) + (m.m[1][2] * p.z) + m.m[1][3];
    T homogeneousZ = (m.m[2][0] * p.x) + (m.m[2][1] * p.y) + (m.m[2][2] * p.z) + m.m[2][3];
    T weight       = (m.m[3][0] * p.x) + (m.m[3][1] * p.y) + (m.m[3][2] * p.z) + m.m[3][3];

   if (weight != 1.0f) {
       return Vector3<T>(
           homogeneousX, 
           homogeneousY,
           homogeneousZ
       ) / weight;
   }
   return Vector3<T>(homogeneousX, homogeneousY, homogeneousZ);
};

template <typename T> inline Vector3<T> Transform::operator()(const Vector3<T> &v) const {
    // No need to use the homogeneous representation of vectors to transform them,
    // because their weight, 0, causes entry (mv)30 to be 0.
    return Vector3<T>(
        (m.m[0][0] * v.x) + (m.m[0][1] * v.y) + (m.m[0][2] * v.z),
        (m.m[1][0] * v.x) + (m.m[1][1] * v.y) + (m.m[1][2] * v.z),
        (m.m[2][0] * v.x) + (m.m[2][1] * v.y) + (m.m[2][2] * v.z)
    );
};

template <typename T> inline Normal3<T> Transform::operator()(const Normal3<T> &n) const {
    // The transform doesn't maintain the orthogonality of the normal and the surface 
    // (the tangent to the surface). The inverse transpose does.
    Matrix4x4 mInvTransp = mInv.Transpose();
    return Normal3<T>(
        (mInvTransp[0][0] * n.x) + (mInvTransp[0][1] * n.y) + (mInvTransp[0][2] * n.z),
        (mInvTransp[1][0] * n.x) + (mInvTransp[1][1] * n.y) + (mInvTransp[1][2] * n.z),
        (mInvTransp[2][0] * n.x) + (mInvTransp[2][1] * n.y) + (mInvTransp[2][2] * n.z)
    );
}

Transform Transform::Translate(const Vector3f &delta) const {
    Matrix4x4 mat(
        1, 0, 0, delta.x,
        0, 1, 0, delta.y,
        0, 0, 1, delta.z,
        0, 0, 0, 1
    );

    Matrix4x4 matInv(
        1, 0, 0, -delta.x,
        0, 1, 0, -delta.y,
        0, 0, 1, -delta.z,
        0, 0, 0, 1
    );

    return Transform(mat, matInv);
}

Transform Transform::Scale(Float sx, Float sy, Float sz) const {
    Matrix4x4 mat(
        sx, 0, 0, 0,
        0, sy, 0, 0,
        0, 0, sz, 0,
        0, 0, 0, 1
    );

    Matrix4x4 matInv(
        1/sx, 0, 0, 0,
        0, 1/sy, 0, 0,
        0, 0, 1/sz, 0,
        0, 0, 0, 1
    );

    return Transform(mat, matInv);
}

Transform Transform::RotateX(Float theta) const {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));

    Matrix4x4 mat(
        1,        0,         0, 0,
        0, cosTheta, -sinTheta, 0,
        0, sinTheta,  cosTheta, 0,
        0,        0,         0, 1
    );

    // A rotation matrix is orthogonal and orthognal matrices have the
    // property that their transpose is also their inverse.
    Matrix4x4 matInv = mat.Transpose();

    return Transform(mat, matInv);
}

Transform Transform::RotateY(Float theta) const {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));

    Matrix4x4 mat(
        cosTheta,  0, sinTheta, 0,
                0,  1,        0, 0,
        -sinTheta,  0, cosTheta, 0,
                0,  0,        0, 1
    );

    // A rotation matrix is orthogonal and orthognal matrices have the
    // property that their transpose is also their inverse.
    Matrix4x4 matInv = mat.Transpose();

    return Transform(mat, matInv);
}

Transform Transform::RotateZ(Float theta) const {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));

    Matrix4x4 mat(
        cosTheta, -sinTheta, 0, 0,
        sinTheta,  cosTheta, 0, 0,
                0,         0, 1, 0,
                0,         0, 0, 1
    );

    // A rotation matrix is orthogonal and orthognal matrices have the
    // property that their transpose is also their inverse.
    Matrix4x4 matInv = mat.Transpose();

    return Transform(mat, matInv);
}

Transform Transform::Rotate(Float theta, const Vector3f & axis) const {
    Float sinTheta = std::sin(theta);
    Float cosTheta = std::cos(theta);

    Vector3f a = Normalize(axis);
    Matrix4x4 mat(
        // Rotation of first standard basis vector.
        a.x * a.x + (1 - a.x * a.x) * cosTheta,
        a.x * a.y * (1 - cosTheta) - a.z * sinTheta,
        a.x * a.z * (1 - cosTheta) + a.y * sinTheta,
        0,

        // Rotation of second standard basis vector.
        a.x * a.y * (1 - cosTheta) + a.z * sinTheta,
        a.y * a.y + (1 - a.y * a.y) * cosTheta,
        a.y * a.z * (1 - cosTheta) - a.x * sinTheta,
        0,

        // Rotation of third standard basis vector.
        a.x * a.z * (1 - cosTheta) - a.y * sinTheta,
        a.y * a.z * (1 - cosTheta) + a.x * sinTheta,
        a.z * a.z + (1 - a.z * a.z) * cosTheta,
        0,

        // Last row.
        0, 0, 0, 1
    );

    return Transform(mat, mat.Transpose());
}

// Transforms a point in camera space to world space. The camera space's frame is
// (pos, ?, up, look), where the iHat basis vector is implicitly defined by the others. 
Transform Transform::LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up) const {
    Matrix4x4 cameraToWorld;

    Vector3f forward = Normalize(look - pos);
    Vector3f left = Normalize(Cross(Normalize(up), forward));
    Vector3f newUp = Cross(forward, left);

    // 1st column: iHat in world coordinates.
    cameraToWorld.m[0][0] = left.x;
    cameraToWorld.m[1][0] = left.y;
    cameraToWorld.m[2][0] = left.z;
    cameraToWorld.m[3][0] = 0.f;

    // 2nd column: jHat in world coordinates.
    cameraToWorld.m[0][1] = newUp.x;
    cameraToWorld.m[1][1] = newUp.y;
    cameraToWorld.m[2][1] = newUp.z;
    cameraToWorld.m[3][1] = 0.f;

    // 3rd column: kHat in world coordinates.
    cameraToWorld.m[0][2] = forward.x;
    cameraToWorld.m[1][2] = forward.y;
    cameraToWorld.m[2][2] = forward.z;
    cameraToWorld.m[3][2] - 0.f;

    // 4th column.
    cameraToWorld.m[0][3] = pos.x;
    cameraToWorld.m[1][3] = pos.y;
    cameraToWorld.m[2][3] = pos.z;
    cameraToWorld.m[3][3] = 1.f;

    return Transform(cameraToWorld.Inverse(), cameraToWorld);
}