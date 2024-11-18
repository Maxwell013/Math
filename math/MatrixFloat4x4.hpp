#pragma once
// MatrixFloat4x4.hpp

#include <iomanip>

#include "VectorFloat3.hpp"

#include "../../Logger/logger.hpp"

namespace math {

    class alignas(16) MatrixFloat4x4 {

    private:

        float m_array[16];

    public:

        static MatrixFloat4x4 Indentity() { return MatrixFloat4x4(); }
        static MatrixFloat4x4 Zero() {
            return MatrixFloat4x4((float[16]) {
                0.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 0.0f
            });
        }

        static MatrixFloat4x4 Scale(const float p_float) {
            return MatrixFloat4x4((float[16]) {
                p_float, 0.0f, 0.0f, 0.0f,
                0.0f, p_float, 0.0f, 0.0f,
                0.0f, 0.0f, p_float, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat4x4 Scale(const float p_x, const float p_y, const float p_z) {
            return MatrixFloat4x4((float[16]) {
                p_x, 0.0f, 0.0f, 0.0f,
                0.0f, p_y, 0.0f, 0.0f,
                0.0f, 0.0f, p_z, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat4x4 Scale(const VectorFloat3& p_vector) {
            return MatrixFloat4x4((float[16]) {
                p_vector.x, 0.0f, 0.0f, 0.0f,
                0.0f, p_vector.y, 0.0f, 0.0f,
                0.0f, 0.0f, p_vector.z, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat4x4 Translation(const float p_x, const float p_y, const float p_z) {
            return MatrixFloat4x4((float[16]) {
                1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                p_x, p_y, p_z, 1.0f
            });
        }

        static MatrixFloat4x4 Translation(const VectorFloat3& p_vector) {
            return MatrixFloat4x4((float[16]) {
                1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                p_vector.x, p_vector.y, p_vector.z, 1.0f
            });
        }

        static MatrixFloat4x4 RotationX(const float p_r) {
            return MatrixFloat4x4((float[16]) {
                1.0f, 0.0f     , 0.0f     , 0.0f,
                0.0f, cosf(p_r),-sinf(p_r), 0.0f,
                0.0f, sinf(p_r), cosf(p_r), 0.0f,
                0.0f, 0.0f     , 0.0f     , 1.0f
            });
        }

        static MatrixFloat4x4 RotationY(const float p_r) {
            return MatrixFloat4x4((float[16]) {
                cosf(p_r), 0.0f, sinf(p_r), 0.0f,
                0.0f     , 1.0f, 0.0f     , 0.0f,
                sinf(p_r), 0.0f, cosf(p_r), 0.0f,
                0.0f     , 0.0f, 0.0f     , 1.0f
            });
        }

        static MatrixFloat4x4 RotationZ(const float p_r) {
            return MatrixFloat4x4((float[16]) {
                cosf(p_r),-sinf(p_r), 0.0f, 0.0f,
                sinf(p_r), cosf(p_r), 0.0f, 0.0f,
                0.0f     , 0.0f     , 1.0f, 0.0f,
                0.0f     , 0.0f     , 0.0f, 1.0f
            });
        }

        static MatrixFloat4x4 Rotation(const float p_r, const VectorFloat3& p_axis) {
            VectorFloat3 unit = p_axis.normalized();
            float x = unit.x;
            float y = unit.y;
            float z = unit.z;
            float s = sinf(p_r);
            float c = cosf(p_r);
            return MatrixFloat4x4((float[16]) {
                (c+x*x*(1-c))  , (x*y*(1-c)-z*s), (x*z*(1-c)-y*s), 0.0f,
                (y*x*(1-c)+z*s), (c+y*y*(1-c))  , (y*z*(1-c)-x*s), 0.0f,
                (z*x*(1-c)-y*s), (z*y*(1-c)+x*s), (c+z*z*(1-c))  , 0.0f,
                0.0f           , 0.0f           , 0.0f           , 1.0f
            });
        }

        static MatrixFloat4x4 Ortho(const float p_left, const float p_right, const float p_bottom, const float p_top, const float p_near, const float p_far) {
            float rl = p_right-p_left;
            float tb = p_top-p_bottom;
            float fn = p_far-p_near;
            return MatrixFloat4x4((float[16]) {
                (2/rl), 0.0f  , 0.0f  , 0.0f,
                0.0f  , (2/tb), 0.0f  , 0.0f,
                0.0f  , 0.0f  ,-(2/fn), 0.0f,
               -(p_right+p_left)/rl, -(p_top+p_bottom)/tb, -(p_far+p_near)/fn, 1.0f
            });
        }

        static MatrixFloat4x4 Perspective(const float p_fov, const float p_aspectRatio, const float p_near, const float p_far) {
            float t = tanf(p_fov/2.0f);
            float fn = p_far-p_near;
            MatrixFloat4x4 result = MatrixFloat4x4((float[16]) {
                (1/p_aspectRatio*t), 0.0f ,  0.0f               ,  0.0f,
                0.0f               , (1/t),  0.0f               ,  0.0f,
                0.0f               , 0.0f , -(p_far+p_near)/fn  , -1.0f,
                0.0f               , 0.0f , -(2*p_far*p_near)/fn,  1.0f
            });
            LOGGER_DEBUG(result);
            return result;
        }

        static MatrixFloat4x4 LookAt(const VectorFloat3& p_position, const VectorFloat3& p_target, const VectorFloat3& p_up) {
            VectorFloat3 direction = (p_position - p_target).normalized();
            VectorFloat3 right = direction.cross(p_up).normalized();
            VectorFloat3 up = right.cross(direction);
            return MatrixFloat4x4((float[16]) {
                right.x    , right.y    , right.z    , 0.0f,
                up.x       , up.y       , up.z       , 0.0f,
                direction.x, direction.y, direction.z, 0.0f,
                0.0f       , 0.0f       , 0.0f       , 1.0f
            }) * Translation(-p_position);
        }

        MatrixFloat4x4() : m_array(
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        ) {}

        MatrixFloat4x4(const float p_array[16]) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] = p_array[i_index];
        }

        MatrixFloat4x4(const MatrixFloat4x4& p_matrix) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] = p_matrix[i_index];
        }

        ~MatrixFloat4x4() {}

        MatrixFloat4x4& operator=(const MatrixFloat4x4& p_matrix) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] = p_matrix[i_index];
            return *this;
        }

        MatrixFloat4x4 operator+(const MatrixFloat4x4& p_matrix) const {
            MatrixFloat4x4 result;
            for (size_t i_index = 0; i_index < 16; i_index++)
                result[i_index] = m_array[i_index] + p_matrix[i_index];
            return result;
        }

        MatrixFloat4x4& operator+=(const MatrixFloat4x4& p_matrix) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] += p_matrix[i_index];
            return *this;
        }

        MatrixFloat4x4 operator-(const MatrixFloat4x4& p_matrix) const {
            MatrixFloat4x4 result;
            for (size_t i_index = 0; i_index < 16; i_index++)
                result[i_index] = m_array[i_index] - p_matrix[i_index];
            return result;
        }

        MatrixFloat4x4& operator-=(const MatrixFloat4x4& p_matrix) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] -= p_matrix[i_index];
            return *this;
        }

        MatrixFloat4x4 operator*(const MatrixFloat4x4& p_matrix) const {
            float array[16];
            for (size_t i_row = 0; i_row < 4; i_row++)
                for (size_t i_column = 0; i_column < 4; i_column++) {
                    array[i_row * 4 + i_column] =
                        m_array[i_column] * p_matrix(i_row, 0) +
                        m_array[4 + i_column] * p_matrix(i_row, 1) +
                        m_array[8 + i_column] * p_matrix(i_row, 2) +
                        m_array[12 + i_column] * p_matrix(i_row, 3);
                }
            return MatrixFloat4x4(array);
        }

        MatrixFloat4x4 operator*(const float p_float) const {
            MatrixFloat4x4 result;
            for (size_t i_index = 0; i_index < 16; i_index++)
                result[i_index] = m_array[i_index] * p_float;
            return result;
        }

        MatrixFloat4x4& operator*=(const float p_float) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] *= p_float;
            return *this;
        }

        MatrixFloat4x4 operator/(const float p_float) const {
            MatrixFloat4x4 result;
            for (size_t i_index = 0; i_index < 16; i_index++)
                result[i_index] = m_array[i_index] / p_float;
            return result;
        }

        MatrixFloat4x4& operator/=(const float p_float) {
            for (size_t i_index = 0; i_index < 16; i_index++)
                m_array[i_index] /= p_float;
            return *this;
        }

        MatrixFloat4x4 transposed() const {
            MatrixFloat4x4 result = *this;
            std::swap<float>(result[1], result[4]);
            std::swap<float>(result[2], result[8]);
            std::swap<float>(result[3], result[12]);
            std::swap<float>(result[6], result[9]);
            std::swap<float>(result[7], result[13]);
            std::swap<float>(result[11], result[14]);
            return result;
        }

        MatrixFloat4x4& transpose() {
            std::swap<float>(m_array[1], m_array[4]);
            std::swap<float>(m_array[2], m_array[8]);
            std::swap<float>(m_array[3], m_array[12]);
            std::swap<float>(m_array[6], m_array[9]);
            std::swap<float>(m_array[7], m_array[13]);
            std::swap<float>(m_array[11], m_array[14]);
            return *this;
        }

        float determinant() const {
            return ( m_array[0] * (
                m_array[5] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[6] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) +
                m_array[7] * (m_array[9] * m_array[14] - m_array[10] * m_array[13])
            ) - m_array[1] * (
                m_array[4] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[6] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[7] * (m_array[8] * m_array[14] - m_array[10] * m_array[12])
            ) + m_array[2] * (
                m_array[4] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) -
                m_array[5] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[7] * (m_array[8] * m_array[13] - m_array[9] * m_array[12])
            ) - m_array[3] * (
                m_array[4] * (m_array[9] * m_array[14] - m_array[10] * m_array[13]) -
                m_array[5] * (m_array[8] * m_array[14] - m_array[10] * m_array[12]) +
                m_array[6] * (m_array[8] * m_array[13] - m_array[9] * m_array[12])
            ));
        }

        MatrixFloat4x4 adjugate() const {
            float array[16] = {
                m_array[5] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[6] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) +
                m_array[7] * (m_array[9] * m_array[14] - m_array[10] * m_array[13]),
                -1.0f * (m_array[4] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[6] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[7] * (m_array[8] * m_array[14] - m_array[10] * m_array[12])),
                m_array[4] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) -
                m_array[5] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[7] * (m_array[8] * m_array[13] - m_array[9] * m_array[12]),
                -1.0f * (m_array[4] * (m_array[9] * m_array[14] - m_array[10] * m_array[13]) -
                m_array[5] * (m_array[8] * m_array[14] - m_array[10] * m_array[12]) +
                m_array[6] * (m_array[8] * m_array[13] - m_array[9] * m_array[12])),
                -1.0f * (m_array[1] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[2] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) +
                m_array[3] * (m_array[9] * m_array[14] - m_array[10] * m_array[13])),
                m_array[0] * (m_array[10] * m_array[15] - m_array[11] * m_array[14]) -
                m_array[2] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[3] * (m_array[8] * m_array[14] - m_array[10] * m_array[12]),
                -1.0f * (m_array[0] * (m_array[9] * m_array[15] - m_array[11] * m_array[13]) -
                m_array[1] * (m_array[8] * m_array[15] - m_array[11] * m_array[12]) +
                m_array[3] * (m_array[8] * m_array[13] - m_array[9] * m_array[12])),
                m_array[0] * (m_array[9] * m_array[14] - m_array[10] * m_array[13]) -
                m_array[1] * (m_array[8] * m_array[14] - m_array[10] * m_array[12]) +
                m_array[2] * (m_array[8] * m_array[13] - m_array[9] * m_array[12]),
                m_array[1] * (m_array[6] * m_array[15] - m_array[7] * m_array[14]) -
                m_array[2] * (m_array[5] * m_array[15] - m_array[7] * m_array[13]) +
                m_array[3] * (m_array[5] * m_array[14] - m_array[6] * m_array[13]),
                -1.0f * (m_array[0] * (m_array[6] * m_array[15] - m_array[7] * m_array[14]) -
                m_array[2] * (m_array[4] * m_array[15] - m_array[7] * m_array[12]) +
                m_array[3] * (m_array[4] * m_array[14] - m_array[6] * m_array[12])),
                m_array[0] * (m_array[5] * m_array[15] - m_array[7] * m_array[13]) -
                m_array[1] * (m_array[4] * m_array[15] - m_array[7] * m_array[12]) +
                m_array[3] * (m_array[4] * m_array[13] - m_array[5] * m_array[12]),
                -1.0f * (m_array[0] * (m_array[5] * m_array[14] - m_array[6] * m_array[13]) -
                m_array[1] * (m_array[4] * m_array[14] - m_array[6] * m_array[12]) +
                m_array[2] * (m_array[4] * m_array[13] - m_array[5] * m_array[12])),
                -1.0f * (m_array[1] * (m_array[6] * m_array[11] - m_array[7] * m_array[10]) -
                m_array[2] * (m_array[5] * m_array[11] - m_array[7] * m_array[9]) +
                m_array[3] * (m_array[5] * m_array[10] - m_array[6] * m_array[9])),
                m_array[0] * (m_array[6] * m_array[11] - m_array[7] * m_array[10]) -
                m_array[2] * (m_array[4] * m_array[11] - m_array[7] * m_array[8]) +
                m_array[3] * (m_array[4] * m_array[10] - m_array[6] * m_array[8]),
                -1.0f * (m_array[0] * (m_array[5] * m_array[11] - m_array[7] * m_array[9]) -
                m_array[1] * (m_array[4] * m_array[11] - m_array[7] * m_array[8]) +
                m_array[3] * (m_array[4] * m_array[9] - m_array[5] * m_array[8])),
                m_array[0] * (m_array[5] * m_array[10] - m_array[6] * m_array[9]) -
                m_array[1] * (m_array[4] * m_array[10] - m_array[6] * m_array[8]) +
                m_array[2] * (m_array[4] * m_array[9] - m_array[5] * m_array[8]),
            };
            return MatrixFloat4x4(array).transpose();
        }

        MatrixFloat4x4 inverse() const {
            const float d = determinant();
            if (abs(d) <= FLT_EPSILON)
                throw std::overflow_error("Divide by zero exception!");
            return adjugate()/d;
        }

        float operator()(const size_t p_row, const size_t p_column) const {
            if (p_row >= 4 || p_column >= 4)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row * 4 + p_column];
        }

        float& operator()(const size_t p_row, const size_t p_column) {
            if (p_row >= 4 || p_column >= 4)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row * 4 + p_column];
        }

        float operator[](const size_t p_index) const {
            if (p_index >= 16)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float& operator[](const size_t p_index) {
            if (p_index >= 16)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float *ptr() {
            return &m_array[0];
        }
    };

    inline bool operator==(const MatrixFloat4x4& p_matrixA, const MatrixFloat4x4& p_matrixB) {
        for (size_t i_index = 0; i_index < 16; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return false;
        return true;
    }

    inline bool operator!=(const MatrixFloat4x4& p_matrixA, const MatrixFloat4x4& p_matrixB) {
        for (size_t i_index = 0; i_index < 16; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return true;
        return false;
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const MatrixFloat4x4& p_matrix) {
        p_ostream << "MatrixFloat4x4 (\n";
        for(size_t i_row = 0; i_row < 4; i_row++)
            p_ostream << "[" << std::fixed << std::setprecision(2)
            << std::setw(5) << p_matrix(i_row, 0) << ", "
            << std::setw(5) << p_matrix(i_row, 1) << ", "
            << std::setw(5) << p_matrix(i_row, 2) << ", "
            << std::setw(5) << p_matrix(i_row, 3) << "]\n";
        p_ostream << ")";
        return p_ostream;
    }
} // namespace math
