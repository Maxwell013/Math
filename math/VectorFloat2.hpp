#pragma once
// VectorFloat2.hpp

#include <iostream>

#include "scalarFloat.hpp"

namespace math {

    class alignas(16) VectorFloat2 {

    public:

        union {
            struct { float x, y; };
            struct { float r, g; };
            struct { float i, j; };
            float components[2];
        };

        static VectorFloat2 Zero() { return VectorFloat2(); };
        static VectorFloat2 X() { return VectorFloat2(1.0f, 0.0f); };
        static VectorFloat2 Y() { return VectorFloat2(0.0f, 1.0f); };

        static VectorFloat2 Right() { return X(); }
        static VectorFloat2 Left() { return Zero() - X(); }
        static VectorFloat2 Up() { return Y(); }
        static VectorFloat2 Down() { return Zero() - Y(); }

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

        VectorFloat2 operator-() const { return VectorFloat2( -x, -y); }

        VectorFloat2& operator-=(const VectorFloat2& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            return *this;
        }

        VectorFloat2 operator*(const VectorFloat2 p_vector) const {
            return VectorFloat2(
                x * p_vector.x,
                y * p_vector.y
            );
        }

        VectorFloat2& operator*=(const VectorFloat2 p_vector) {
            x *= p_vector.x;
            y *= p_vector.y;
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

        VectorFloat2 normalized() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return VectorFloat2(0, 0);
            return VectorFloat2(
                x / m,
                y / m
            );
        }

        VectorFloat2& normalize() {
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
