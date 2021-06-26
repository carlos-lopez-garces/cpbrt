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
    Matrix4x4 mInvTransp = mInv.Transpose();
    return Normal3<T>(
        (mInvTransp[0][0] * n.x) + (mInvTransp[0][1] * n.y) + (mInvTransp[0][2] * n.z),
        (mInvTransp[1][0] * n.x) + (mInvTransp[1][1] * n.y) + (mInvTransp[1][2] * n.z),
        (mInvTransp[2][0] * n.x) + (mInvTransp[2][1] * n.y) + (mInvTransp[2][2] * n.z)
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

SurfaceInteraction Transform::operator()(const SurfaceInteraction &si) const {
    SurfaceInteraction transformedSI;
    // TODO: implement when ready to handle floating-point errors.
    return transformedSI;
}

inline Transform Transform::operator*(const Transform &t) const {
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