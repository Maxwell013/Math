#pragma once
// MatrixFloat2x2.hpp

#include <cfloat>
#include <iomanip>

namespace math {

    class alignas(16) MatrixFloat2x2 {

    private:

        float m_array[4];

    public:

        static MatrixFloat2x2 Indentity() { return MatrixFloat2x2(); }
        static MatrixFloat2x2 Zero() { return MatrixFloat2x2((float[4]) {
            0.0f, 0.0f,
            0.0f, 0.0f
        }); }

        MatrixFloat2x2() : m_array(
            1.0f, 0.0f,
            0.0f, 1.0f
        ) {}

        MatrixFloat2x2(const float p_array[4]) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] = p_array[i_index];
        }

        MatrixFloat2x2(const MatrixFloat2x2 &p_matrix) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] = p_matrix[i_index];
        }

        ~MatrixFloat2x2() {}

        MatrixFloat2x2& operator=(const MatrixFloat2x2& p_matrix) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] = p_matrix[i_index];
            return *this;
        }

        MatrixFloat2x2 operator+(const MatrixFloat2x2& p_matrix) const {
            MatrixFloat2x2 result;
            for (size_t i_index = 0; i_index < 4; i_index++)
                result[i_index] = m_array[i_index] + p_matrix[i_index];
            return result;
        }

        MatrixFloat2x2& operator+=(const MatrixFloat2x2& p_matrix) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] += p_matrix[i_index];
            return *this;
        }

        MatrixFloat2x2 operator-(const MatrixFloat2x2& p_matrix) const {
            MatrixFloat2x2 result;
            for (size_t i_index = 0; i_index < 4; i_index++)
                result[i_index] = m_array[i_index] - p_matrix[i_index];
            return result;
        }

        MatrixFloat2x2& operator-=(const MatrixFloat2x2& p_matrix) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] -= p_matrix[i_index];
            return *this;
        }

        MatrixFloat2x2 operator*(const MatrixFloat2x2& p_matrix) const {
            float array[4];
            for (size_t i_row = 0; i_row < 2; i_row++)
                for (size_t i_column = 0; i_column < 2; i_column++) {
                    array[i_row + 2*i_column] =
                        m_array[i_row] * p_matrix[2*i_column] +
                        m_array[i_row + 2] * p_matrix[1 + 2*i_column];
                }
            return MatrixFloat2x2(array);
        }

        MatrixFloat2x2& operator*=(const MatrixFloat2x2& p_matrix) {
            for (size_t i_row = 0; i_row < 2; i_row++)
                for (size_t i_column = 0; i_column < 2; i_column++) {
                    m_array[i_row + 2*i_column] =
                        m_array[i_row] * p_matrix[2*i_column] +
                        m_array[i_row + 2] * p_matrix[1 + 2*i_column];
                }
            return *this;
        }

        MatrixFloat2x2 operator*(const float p_float) const {
            MatrixFloat2x2 result;
            for (size_t i_index = 0; i_index < 4; i_index++)
                result[i_index] = m_array[i_index] * p_float;
            return result;
        }

        MatrixFloat2x2& operator*=(const float p_float) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] *= p_float;
            return *this;
        }

        MatrixFloat2x2 operator/(const float p_float) const {
            MatrixFloat2x2 result;
            for (size_t i_index = 0; i_index < 4; i_index++)
                result[i_index] = m_array[i_index] / p_float;
            return result;
        }

        MatrixFloat2x2& operator/=(const float p_float) {
            for (size_t i_index = 0; i_index < 4; i_index++)
                m_array[i_index] /= p_float;
            return *this;
        }

        MatrixFloat2x2 transposed() const {
            MatrixFloat2x2 result = *this;
            std::swap<float>(result[1], result[2]);
            return result;
        }

        MatrixFloat2x2& transpose() {
            std::swap<float>(m_array[1], m_array[2]);
            return *this;
        }

        float determinant() const {
            return m_array[0] * m_array[3] - m_array[1] * m_array[2];
        }

        MatrixFloat2x2 adjugate() const {
            float array[4] = {m_array[3], m_array[2], m_array[1], m_array[0]};
            return MatrixFloat2x2(array).transpose();
        }

        MatrixFloat2x2 inverse() const {
            const float d = determinant();
            if (abs(d) <= FLT_EPSILON)
                throw std::overflow_error("Divide by zero exception!");
            return adjugate()/d;
        }

        float operator()(const size_t p_row, const size_t p_column) const {
            if (p_row >= 2 || p_column >= 2)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row + p_column * 2];
        }

        float& operator()(const size_t p_row, const size_t p_column) {
            if (p_row >= 2 || p_column >= 2)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_row + p_column * 2];
        }

        float operator[](const size_t p_index) const {
            if (p_index >= 4)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float& operator[](const size_t p_index) {
            if (p_index >= 4)
                throw std::out_of_range("Matrix index out of range!");
            return m_array[p_index];
        }

        float *ptr() {
            return &m_array[0];
        }
    };

    inline bool operator==(const MatrixFloat2x2& p_matrixA, const MatrixFloat2x2& p_matrixB) {
        for (size_t i_index = 0; i_index < 4; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return false;
        return true;
    }

    inline bool operator!=(const MatrixFloat2x2& p_matrixA, const MatrixFloat2x2& p_matrixB) {
        for (size_t i_index = 0; i_index < 4; i_index++)
            if (abs(p_matrixA[i_index] - p_matrixB[i_index]) > FLT_EPSILON) return true;
        return false;
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const MatrixFloat2x2& p_matrix) {
        p_ostream << "MatrixFloat2x2 (\n";
        for(size_t i_row = 0; i_row < 2; i_row++)
            p_ostream << "[" << std::fixed << std::setprecision(2)
            << std::setw(5) << p_matrix(i_row, 0) << ", "
            << std::setw(5) << p_matrix(i_row, 1) << "]\n";
        p_ostream << ")";
        return p_ostream;
    }
} // namespace math
