#pragma once
// MatrixFloat4x4.hpp

#include <cfloat>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cstddef>
#include <stdexcept>

#include "../simd.hpp"

namespace math {

    class alignas(16) MatrixFloat4x4 {

    private:

        float m_array[16];

    public:

        static MatrixFloat4x4 Indentity() { return MatrixFloat4x4(); }
        static MatrixFloat4x4 Zero() { return MatrixFloat4x4((float[16]) {
            0.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 0.0f
        }); }

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

        MatrixFloat4x4(const MatrixFloat4x4 &p_matrix) {
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

            SIMD_F32X4_T A0;
            SIMD_F32X4_T A1;
            SIMD_F32X4_T A2;
            SIMD_F32X4_T A3;

            SIMD_F32X4_T B0;
            SIMD_F32X4_T B1;
            SIMD_F32X4_T B2;
            SIMD_F32X4_T B3;

            SIMD_F32X4_T C0;
            SIMD_F32X4_T C1;
            SIMD_F32X4_T C2;
            SIMD_F32X4_T C3;

            A0 = SIMD_LOAD_F32(m_array);
            A1 = SIMD_LOAD_F32(m_array + 4);
            A2 = SIMD_LOAD_F32(m_array + 8);
            A3 = SIMD_LOAD_F32(m_array + 12);

            C0 = SIMD_MOVE_F32(0);
            C1 = SIMD_MOVE_F32(0);
            C2 = SIMD_MOVE_F32(0);
            C3 = SIMD_MOVE_F32(0);

            B0 = SIMD_LOAD_F32(p_matrix.m_array);
            C0 = SIMD_FMA_F32(C0, A0, B0, 0);
            C0 = SIMD_FMA_F32(C0, A1, B0, 1);
            C0 = SIMD_FMA_F32(C0, A2, B0, 2);
            C0 = SIMD_FMA_F32(C0, A3, B0, 3);
            SIMD_STORE_F32(array, C0);

            B1 = SIMD_LOAD_F32(p_matrix.m_array + 4);
            C0 = SIMD_FMA_F32(C0, A0, B0, 0);
            C1 = SIMD_FMA_F32(C1, A0, B1, 0);
            C1 = SIMD_FMA_F32(C1, A1, B1, 1);
            C1 = SIMD_FMA_F32(C1, A2, B1, 2);
            C1 = SIMD_FMA_F32(C1, A3, B1, 3);
            SIMD_STORE_F32(array + 4, C1);

            B2 = SIMD_LOAD_F32(p_matrix.m_array + 8);
            C2 = SIMD_FMA_F32(C2, A0, B2, 0);
            C2 = SIMD_FMA_F32(C2, A1, B2, 1);
            C2 = SIMD_FMA_F32(C2, A2, B2, 2);
            C2 = SIMD_FMA_F32(C2, A3, B2, 3);
            SIMD_STORE_F32(array + 8, C2);

            B3 = SIMD_LOAD_F32(p_matrix.m_array + 12);
            C3 = SIMD_FMA_F32(C3, A0, B3, 0);
            C3 = SIMD_FMA_F32(C3, A1, B3, 1);
            C3 = SIMD_FMA_F32(C3, A2, B3, 2);
            C3 = SIMD_FMA_F32(C3, A3, B3, 3);
            SIMD_STORE_F32(array + 12, C3);

            return MatrixFloat4x4(array);
        }

        MatrixFloat4x4& operator*=(const MatrixFloat4x4& p_matrix) {
            SIMD_F32X4_T A0;
            SIMD_F32X4_T A1;
            SIMD_F32X4_T A2;
            SIMD_F32X4_T A3;

            SIMD_F32X4_T B0;
            SIMD_F32X4_T B1;
            SIMD_F32X4_T B2;
            SIMD_F32X4_T B3;

            SIMD_F32X4_T C0;
            SIMD_F32X4_T C1;
            SIMD_F32X4_T C2;
            SIMD_F32X4_T C3;

            A0 = SIMD_LOAD_F32(m_array);
            A1 = SIMD_LOAD_F32(m_array + 4);
            A2 = SIMD_LOAD_F32(m_array + 8);
            A3 = SIMD_LOAD_F32(m_array + 12);

            C0 = SIMD_MOVE_F32(0);
            C1 = SIMD_MOVE_F32(0);
            C2 = SIMD_MOVE_F32(0);
            C3 = SIMD_MOVE_F32(0);

            B0 = SIMD_LOAD_F32(p_matrix.m_array);
            C0 = SIMD_FMA_F32(C0, A0, B0, 0);
            C0 = SIMD_FMA_F32(C0, A1, B0, 1);
            C0 = SIMD_FMA_F32(C0, A2, B0, 2);
            C0 = SIMD_FMA_F32(C0, A3, B0, 3);
            SIMD_STORE_F32(m_array, C0);

            B1 = SIMD_LOAD_F32(p_matrix.m_array + 4);
            C0 = SIMD_FMA_F32(C0, A0, B0, 0);
            C1 = SIMD_FMA_F32(C1, A0, B1, 0);
            C1 = SIMD_FMA_F32(C1, A1, B1, 1);
            C1 = SIMD_FMA_F32(C1, A2, B1, 2);
            C1 = SIMD_FMA_F32(C1, A3, B1, 3);
            SIMD_STORE_F32(m_array + 4, C1);

            B2 = SIMD_LOAD_F32(p_matrix.m_array + 8);
            C2 = SIMD_FMA_F32(C2, A0, B2, 0);
            C2 = SIMD_FMA_F32(C2, A1, B2, 1);
            C2 = SIMD_FMA_F32(C2, A2, B2, 2);
            C2 = SIMD_FMA_F32(C2, A3, B2, 3);
            SIMD_STORE_F32(m_array + 8, C2);

            B3 = SIMD_LOAD_F32(p_matrix.m_array + 12);
            C3 = SIMD_FMA_F32(C3, A0, B3, 0);
            C3 = SIMD_FMA_F32(C3, A1, B3, 1);
            C3 = SIMD_FMA_F32(C3, A2, B3, 2);
            C3 = SIMD_FMA_F32(C3, A3, B3, 3);
            SIMD_STORE_F32(m_array + 12, C3);

            return *this;
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
            return m_array[p_row + p_column * 4];
        }

        float& operator()(const size_t p_row, const size_t p_column) {
            if (p_row >= 4 || p_column >= 4)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row + p_column * 4];
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
