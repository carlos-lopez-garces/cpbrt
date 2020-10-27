#include "transform.h"

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