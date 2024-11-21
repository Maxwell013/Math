#pragma once
// MatrixFloat3x3.hpp

#include <iomanip>

#include "VectorFloat2.hpp"
#include "trigonometryFloat.hpp"

namespace math {

    class alignas(16) MatrixFloat3x3 {

    private:

        float m_array[9];

    public:

        static MatrixFloat3x3 Indentity() { return MatrixFloat3x3(); }
        static MatrixFloat3x3 Zero() { return MatrixFloat3x3((float[9]) {
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f
        }); }

        static MatrixFloat3x3 Scale(const float p_float) {
            return MatrixFloat3x3((float[9]) {
                p_float, 0.0f, 0.0f,
                0.0f, p_float, 0.0f,
                0.0f, 0.0f, 1.0f,
            });
        }

        static MatrixFloat3x3 Scale(const float p_x, const float p_y) {
            return MatrixFloat3x3((float[9]) {
                p_x, 0.0f, 0.0f,
                0.0f, p_y, 0.0f,
                0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat3x3 Scale(const VectorFloat2& p_vector) {
            return MatrixFloat3x3((float[9]) {
                p_vector.x, 0.0f, 0.0f,
                0.0f, p_vector.y, 0.0f,
                0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat3x3 Translation(const float p_x, const float p_y) {
            return MatrixFloat3x3((float[9]) {
                1.0f, 0.0f, p_x,
                0.0f, 1.0f, p_y,
                0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat3x3 Translation(const VectorFloat2& p_vector) {
            return MatrixFloat3x3((float[9]) {
                1.0f, 0.0f, p_vector.x,
                0.0f, 1.0f, p_vector.y,
                0.0f, 0.0f, 1.0f
            });
        }

        static MatrixFloat3x3 Rotation(const float p_r) {
                return MatrixFloat3x3((float[9]) {
                cos(p_r),-sin(p_r), 0.0f,
                sin(p_r), cos(p_r), 0.0f,
                0.0f    , 0.0f    , 1.0f
            });
        }

        MatrixFloat3x3() : m_array(
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        ) {}

        MatrixFloat3x3(const float p_array[9]) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] = p_array[i_index];
        }

        MatrixFloat3x3(const MatrixFloat3x3 &p_matrix) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] = p_matrix[i_index];
        }

        ~MatrixFloat3x3() {}

        MatrixFloat3x3& operator=(const MatrixFloat3x3& p_matrix) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] = p_matrix[i_index];
            return *this;
        }

        MatrixFloat3x3 operator+(const MatrixFloat3x3& p_matrix) const {
            MatrixFloat3x3 result;
            for (size_t i_index = 0; i_index < 9; i_index++)
                result[i_index] = m_array[i_index] + p_matrix[i_index];
            return result;
        }

        MatrixFloat3x3& operator+=(const MatrixFloat3x3& p_matrix) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] += p_matrix[i_index];
            return *this;
        }

        MatrixFloat3x3 operator-(const MatrixFloat3x3& p_matrix) const {
            MatrixFloat3x3 result;
            for (size_t i_index = 0; i_index < 9; i_index++)
                result[i_index] = m_array[i_index] - p_matrix[i_index];
            return result;
        }

        MatrixFloat3x3& operator-=(const MatrixFloat3x3& p_matrix) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] -= p_matrix[i_index];
            return *this;
        }

        MatrixFloat3x3 operator*(const MatrixFloat3x3& p_matrix) const {
            float array[9];
            for (size_t i_row = 0; i_row < 3; i_row++)
                for (size_t i_column = 0; i_column < 3; i_column++) {
                    array[i_row + 3*i_column] =
                        m_array[i_row] * p_matrix[3*i_column] +
                        m_array[i_row + 3] * p_matrix[1 + 3*i_column] +
                        m_array[i_row + 6] * p_matrix[2 + 3*i_column];
                }
            return MatrixFloat3x3(array);
        }

        MatrixFloat3x3 operator*(const float p_float) const {
            MatrixFloat3x3 result;
            for (size_t i_index = 0; i_index < 9; i_index++)
                result[i_index] = m_array[i_index] * p_float;
            return result;
        }

        MatrixFloat3x3& operator*=(const float p_float) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] *= p_float;
            return *this;
        }

        MatrixFloat3x3 operator/(const float p_float) const {
            MatrixFloat3x3 result;
            for (size_t i_index = 0; i_index < 9; i_index++)
                result[i_index] = m_array[i_index] / p_float;
            return result;
        }

        MatrixFloat3x3& operator/=(const float p_float) {
            for (size_t i_index = 0; i_index < 9; i_index++)
                m_array[i_index] /= p_float;
            return *this;
        }

        MatrixFloat3x3 transposed() const {
            MatrixFloat3x3 result = *this;
            std::swap<float>(result[1], result[3]);
            std::swap<float>(result[2], result[6]);
            std::swap<float>(result[5], result[7]);
            return result;
        }

        MatrixFloat3x3& transpose() {
            std::swap<float>(m_array[1], m_array[3]);
            std::swap<float>(m_array[2], m_array[6]);
            std::swap<float>(m_array[5], m_array[7]);
            return *this;
        }

        float determinant() const {
            return (
                m_array[0] * (m_array[4] * m_array[8] - m_array[5] * m_array[7]) -
                m_array[1] * (m_array[3] * m_array[8] - m_array[5] * m_array[6]) +
                m_array[2] * (m_array[3] * m_array[7] - m_array[4] * m_array[6])
            );
        }

        MatrixFloat3x3 adjugate() const {
            float array[9] = {
                m_array[4] * m_array[8] - m_array[5] * m_array[7],
                m_array[3] * m_array[8] - m_array[5] * m_array[6],
                m_array[3] * m_array[7] - m_array[4] * m_array[6],
                m_array[1] * m_array[8] - m_array[2] * m_array[7],
                m_array[0] * m_array[8] - m_array[2] * m_array[6],
                m_array[0] * m_array[7] - m_array[1] * m_array[6],
                m_array[1] * m_array[5] - m_array[2] * m_array[4],
                m_array[0] * m_array[5] - m_array[3] * m_array[2],
                m_array[0] * m_array[4] - m_array[1] * m_array[3]
            };
            return MatrixFloat3x3(array).transpose();
        }

        MatrixFloat3x3 inverse() const {
            const float d = determinant();
            if (abs(d) <= FLT_EPSILON)
                throw std::overflow_error("Divide by zero exception!");
            return adjugate()/d;
        }

        float operator()(const size_t p_row, const size_t p_column) const {
            if (p_row >= 3 || p_column >= 3)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row * 3 + p_column];
        }

        float& operator()(const size_t p_row, const size_t p_column) {
            if (p_row >= 3 || p_column >= 3)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row * 3 + p_column];
        }

        float operator[](const size_t p_index) const {
            if (p_index >= 9)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float& operator[](const size_t p_index) {
            if (p_index >= 9)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float *ptr() {
            return &m_array[0];
        }
    };

    inline bool operator==(const MatrixFloat3x3& p_matrixA, const MatrixFloat3x3& p_matrixB) {
        for (size_t i_index = 0; i_index < 9; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return false;
        return true;
    }

    inline bool operator!=(const MatrixFloat3x3& p_matrixA, const MatrixFloat3x3& p_matrixB) {
        for (size_t i_index = 0; i_index < 9; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return true;
        return false;
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const MatrixFloat3x3& p_matrix) {
        p_ostream << "MatrixFloat3x3 (\n";
        for(size_t i_row = 0; i_row < 3; i_row++)
            p_ostream << "[" << std::fixed << std::setprecision(2)
            << std::setw(5) << p_matrix(i_row, 0) << ", "
            << std::setw(5) << p_matrix(i_row, 1) << ", "
            << std::setw(5) << p_matrix(i_row, 2) << "]\n";
        p_ostream << ")";
        return p_ostream;
    }
} // namespace math
