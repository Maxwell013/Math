#pragma once
// VectorFloat3.hpp

#include <iostream>

#include "scalarFloat.hpp"

namespace math {

    class alignas(16) VectorFloat3 {

    public:

        union {
            struct { float x, y, z; };
            struct { float r, g, b; };
            struct { float i, j, k; };
            float components[3];
        };

        static VectorFloat3 Zero() { return VectorFloat3(0.0f, 0.0f, 0.0f); }
        static VectorFloat3 X() { return VectorFloat3(1.0f, 0.0f, 0.0f); }
        static VectorFloat3 Y() { return VectorFloat3(0.0f, 1.0f, 0.0f); }
        static VectorFloat3 Z() { return VectorFloat3(0.0f, 0.0f, 1.0f); }

        static VectorFloat3 Right() { return X(); }
        static VectorFloat3 Left() { return Zero() - X(); }
        static VectorFloat3 Up() { return Y(); }
        static VectorFloat3 Down() { return Zero() - Y(); }
        static VectorFloat3 Foward() { return Z(); }
        static VectorFloat3 Backward() { return Zero() - Z(); }

        VectorFloat3(const float p_x = 0.0f, const float p_y = 0.0f, const float p_z = 1.0f) : x(p_x), y(p_y), z(p_z) {}

        VectorFloat3(const VectorFloat3& p_vector) : x(p_vector.x), y(p_vector.y), z(p_vector.z) {}

        ~VectorFloat3() {}

        VectorFloat3& operator=(const VectorFloat3& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            z = p_vector.z;
            return *this;
        }

        VectorFloat3 operator+(const VectorFloat3& p_vector) const {
            return VectorFloat3(
                x + p_vector.x,
                y + p_vector.y,
                z + p_vector.z
            );
        }

        VectorFloat3& operator+=(const VectorFloat3& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            z += p_vector.z;
            return *this;
        }

        VectorFloat3 operator-(const VectorFloat3& p_vector) const {
            return VectorFloat3(
                x - p_vector.x,
                y - p_vector.y,
                z - p_vector.z
            );
        }

        VectorFloat3 operator-() const { return VectorFloat3( -x, -y, -z); }

        VectorFloat3& operator-=(const VectorFloat3& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            z -= p_vector.z;
            return *this;
        }

        VectorFloat3 operator*(const VectorFloat3 p_vector) const {
            return VectorFloat3(
                x * p_vector.x,
                y * p_vector.y,
                z * p_vector.z
            );
        }

        VectorFloat3& operator*=(const VectorFloat3 p_vector) {
            x *= p_vector.x;
            y *= p_vector.y;
            z *= p_vector.z;
            return *this;
        }

        VectorFloat3 operator*(const float p_float) const {
            return VectorFloat3(
                x * p_float,
                y * p_float,
                z * p_float
            );
        }

        VectorFloat3& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            z *= p_float;
            return *this;
        }

        VectorFloat3 operator/(const float p_float) const {
            return VectorFloat3(
                x / p_float,
                y / p_float,
                z / p_float
            );
        }

        VectorFloat3& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            z /= p_float;
            return *this;
        }

        float dot(const VectorFloat3& p_vector) const {
            return (
            x * p_vector.x +
            y * p_vector.y +
            z * p_vector.z
            );
        }

        VectorFloat3 cross(const VectorFloat3& p_vector) const {
            return VectorFloat3(
                y * p_vector.z - z * p_vector.y,
                z * p_vector.x - x * p_vector.z,
                x * p_vector.y - y * p_vector.x
            );
        }

        float magnitude() const {
            return sqrt(
                x*x +
                y*y +
                z*z
            );
        }

        VectorFloat3 normalized() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return VectorFloat3(0, 0, 0);
            return VectorFloat3(
                x / m,
                y / m,
                z / m
            );
        }

        VectorFloat3& normalize() {
            float m = magnitude();
            if (m <= FLT_EPSILON) {
                x = 0; y = 0; z = 0;
            } else {
                x /= m;
                y /= m;
                z /= m;
            }
            return *this;
        }
    };

    inline bool operator==(const VectorFloat3& p_vectorA, const VectorFloat3& p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) <= FLT_EPSILON &&
            abs(p_vectorA.y - p_vectorB.y) <= FLT_EPSILON &&
            abs(p_vectorA.z - p_vectorB.z) <= FLT_EPSILON
        );
    }

    inline bool operator!=(const VectorFloat3 &p_vectorA, const VectorFloat3 &p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) > FLT_EPSILON ||
            abs(p_vectorA.y - p_vectorB.y) > FLT_EPSILON ||
            abs(p_vectorA.z - p_vectorB.z) > FLT_EPSILON
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const VectorFloat3& p_vector) {
        return p_ostream << "VectorFloat3(" << p_vector.x << ", " << p_vector.y << ", " << p_vector.z << ")";
    }

} // namespace math
