#include "transform.h"
#include "interaction.h"

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

Matrix4x4 Matrix4x4::Inverse() const {
    return ::Inverse(*this);
}

Matrix4x4 Inverse(const Matrix4x4 &m) {
    int indxc[4], indxr[4];
    int ipiv[4] = {0, 0, 0, 0};
    Float minv[4][4];
    memcpy(minv, m.m, 4 * 4 * sizeof(Float));
    for (int i = 0; i < 4; i++) {
        int irow = 0, icol = 0;
        Float big = 0.f;
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(minv[j][k]) >= big) {
                            big = Float(std::abs(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1)
                        Error("Singular matrix in MatrixInvert");
                }
            }
        }
        ++ipiv[icol];
        if (irow != icol) {
            for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (minv[icol][icol] == 0.f) Error("Singular matrix in MatrixInvert");

        Float pivinv = 1. / minv[icol][icol];
        minv[icol][icol] = 1.;
        for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                Float save = minv[j][icol];
                minv[j][icol] = 0;
                for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
            }
        }
    }

    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
        }
    }

    return Matrix4x4(minv);
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

SurfaceInteraction Transform::operator()(const SurfaceInteraction &si) const {
    SurfaceInteraction ret;
    ret.p = (*this)(si.p, si.pError, &ret.pError);

    const Transform &t = *this;
    ret.n = Normalize(t(si.n));
    ret.wo = Normalize(t(si.wo));
    ret.time = si.time;
    ret.mediumInterface = si.mediumInterface;
    ret.uv = si.uv;
    ret.shape = si.shape;
    ret.dpdu = t(si.dpdu);
    ret.dpdv = t(si.dpdv);
    ret.dndu = t(si.dndu);
    ret.dndv = t(si.dndv);
    ret.shading.n = Normalize(t(si.shading.n));
    ret.shading.dpdu = t(si.shading.dpdu);
    ret.shading.dpdv = t(si.shading.dpdv);
    ret.shading.dndu = t(si.shading.dndu);
    ret.shading.dndv = t(si.shading.dndv);
    ret.dudx = si.dudx;
    ret.dvdx = si.dvdx;
    ret.dudy = si.dudy;
    ret.dvdy = si.dvdy;
    ret.dpdx = t(si.dpdx);
    ret.dpdy = t(si.dpdy);
    ret.bsdf = si.bsdf;
    ret.bssrdf = si.bssrdf;
    ret.primitive = si.primitive;
    ret.shading.n = Faceforward(ret.shading.n, ret.n);
    return ret;
}

Transform Transform::operator*(const Transform &t) const {
    // The inverse of the product is the product of the inverses in reverse order.
    return Transform(m*t.m, t.mInv*mInv);
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
    cameraToWorld.m[3][2] = 0.f;

    // 4th column.
    cameraToWorld.m[0][3] = pos.x;
    cameraToWorld.m[1][3] = pos.y;
    cameraToWorld.m[2][3] = pos.z;
    cameraToWorld.m[3][3] = 1.f;

    return Transform(cameraToWorld.Inverse(), cameraToWorld);
}

bool Transform::SwapsHandedness() const {
    // Compute the determinant of the upper-left 3x3 submatrix using a cofactor
    // expansion across the first row.
    Float upperLeft3x3Determinant = 
        m.m[0][0]   * (m.m[1][1]*m.m[2][2] - m.m[1][2]*m.m[2][1])
        - m.m[0][1] * (m.m[1][0]*m.m[2][2] - m.m[1][2]*m.m[2][0])
        + m.m[0][2] * (m.m[1][0]*m.m[2][1] - m.m[1][1]*m.m[2][0]);

    // The transformation changes the handedness of the coordinate system when this
    // determinant is negative.
    return upperLeft3x3Determinant < 0.f;
}

Transform Translate(const Vector3f &delta) {
    Matrix4x4 m(1, 0, 0, delta.x, 0, 1, 0, delta.y, 0, 0, 1, delta.z, 0, 0, 0, 1);
    Matrix4x4 minv(1, 0, 0, -delta.x, 0, 1, 0, -delta.y, 0, 0, 1, -delta.z, 0, 0, 0, 1);
    return Transform(m, minv);
}

Transform Scale(Float x, Float y, Float z) {
    Matrix4x4 m(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
    Matrix4x4 minv(1 / x, 0, 0, 0, 0, 1 / y, 0, 0, 0, 0, 1 / z, 0, 0, 0, 0, 1);
    return Transform(m, minv);
}

Transform RotateX(Float theta) {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(1, 0, 0, 0, 0, cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0, 0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateZ(Float theta) {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform Rotate(Float theta, const Vector3f &axis) {
    Vector3f a = Normalize(axis);
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    Matrix4x4 m;

    // Compute rotation of first basis vector.
    m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
    m.m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
    m.m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
    m.m[0][3] = 0;

    // Compute rotations of second and third basis vectors.
    m.m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
    m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
    m.m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
    m.m[1][3] = 0;

    m.m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
    m.m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
    m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
    m.m[2][3] = 0;

    return Transform(m, Transpose(m));
}

Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up) {
    Matrix4x4 cameraToWorld;

    // Initialize fourth column of viewing matrix.
    cameraToWorld.m[0][3] = pos.x;
    cameraToWorld.m[1][3] = pos.y;
    cameraToWorld.m[2][3] = pos.z;
    cameraToWorld.m[3][3] = 1;

    // Initialize first three columns of viewing matrix.
    Vector3f dir = Normalize(look - pos);
    if (Cross(Normalize(up), dir).Length() == 0) {
        Error(
            "\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) "
            "passed to LookAt are pointing in the same direction.  Using "
            "the identity transformation.",
            up.x, up.y, up.z, dir.x, dir.y, dir.z);
        return Transform();
    }

    Vector3f right = Normalize(Cross(Normalize(up), dir));

    Vector3f newUp = Cross(dir, right);

    cameraToWorld.m[0][0] = right.x;
    cameraToWorld.m[1][0] = right.y;
    cameraToWorld.m[2][0] = right.z;
    cameraToWorld.m[3][0] = 0.;
    cameraToWorld.m[0][1] = newUp.x;
    cameraToWorld.m[1][1] = newUp.y;
    cameraToWorld.m[2][1] = newUp.z;
    cameraToWorld.m[3][1] = 0.;
    cameraToWorld.m[0][2] = dir.x;
    cameraToWorld.m[1][2] = dir.y;
    cameraToWorld.m[2][2] = dir.z;
    cameraToWorld.m[3][2] = 0.;

    return Transform(Inverse(cameraToWorld), cameraToWorld);
}

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