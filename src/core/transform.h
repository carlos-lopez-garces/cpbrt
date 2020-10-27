#include "pbrt.h"

class Matrix4x4 {
private:
    Float m[4][4];

public:
    Matrix4x4() {
        // Identity matrix.
        m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
        m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
        m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
        m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
    }

    Matrix4x4(Float m[4][4]) {
        memcpy(this->m, m, 16 * sizeof(Float));
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
};