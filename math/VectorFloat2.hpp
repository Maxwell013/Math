// VectorFloat2.hpp
#pragma once

#include <cfloat>
#include <cmath>
#include <iostream>

namespace math {

    class VectorFloat2 {

    public:

        float x, y;

        VectorFloat2(const float p_x = 0, const float p_y = 0) : x(p_x), y(p_y) {}

        VectorFloat2(const VectorFloat2& p_vector) : x(p_vector.x), y(p_vector.y) {}

        ~VectorFloat2() {}

        VectorFloat2& operator=(const VectorFloat2& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            return *this;
        }

        VectorFloat2 operator+(const VectorFloat2& p_vector) const {
            return VectorFloat2(
                x + p_vector.x,
                y + p_vector.y
            );
        }

        VectorFloat2& operator+=(const VectorFloat2& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            return *this;
        }

        VectorFloat2 operator-(const VectorFloat2& p_vector) const {
            return VectorFloat2(
                x - p_vector.x,
                y - p_vector.y
            );
        }

        VectorFloat2& operator-=(const VectorFloat2& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            return *this;
        }

        VectorFloat2 operator*(const float p_float) const {
            return VectorFloat2(
                x * p_float,
                y * p_float
            );
        }

        VectorFloat2& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            return *this;
        }

        VectorFloat2 operator/(const float p_float) const {
            return VectorFloat2(
                x / p_float,
                y / p_float
            );
        }

        VectorFloat2& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            return *this;
        }

        float dot(const VectorFloat2& p_vector) const {
            return (
            x * p_vector.x +
            y * p_vector.y
            );
        }

        float magnitude() const {
            return sqrt(
                x*x +
                y*y
            );
        }

        VectorFloat2 normalize() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return VectorFloat2(0, 0);
            return VectorFloat2(
                x / m,
                y / m
            );
        }

        VectorFloat2& normalized() {
            float m = magnitude();
            if (m <= FLT_EPSILON) {
                x = 0; y = 0;
            } else {
                x /= m;
                y /= m;
            }
            return *this;
        }
    };

    inline bool operator==(const VectorFloat2& p_vectorA, const VectorFloat2& p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) <= FLT_EPSILON &&
            abs(p_vectorA.y - p_vectorB.y) <= FLT_EPSILON
        );
    }

    inline bool operator!=(const VectorFloat2 &p_vectorA, const VectorFloat2 &p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) > FLT_EPSILON ||
            abs(p_vectorA.y - p_vectorB.y) > FLT_EPSILON
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const VectorFloat2& p_vector) {
        return p_ostream << "VectorFloat2(" << p_vector.x << ", " << p_vector.y << ")";
    }

} // namespace math
