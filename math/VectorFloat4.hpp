// VectorFloat4.hpp
#pragma once

#include <cfloat>
#include <cmath>
#include <iostream>

namespace math {

    class VectorFloat4 {

    public:

        float x, y, z, w;

        static VectorFloat4 Zero() { return VectorFloat4(); }

        VectorFloat4(const float p_x = 0, const float p_y = 0, const float p_z = 0, const float p_w = 0) : x(p_x), y(p_y), z(p_z), w(p_w) {}

        VectorFloat4(const VectorFloat4& p_vector) : x(p_vector.x), y(p_vector.y), z(p_vector.z), w(p_vector.w) {}

        ~VectorFloat4() {}

        VectorFloat4& operator=(const VectorFloat4& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            z = p_vector.z;
            w = p_vector.w;
            return *this;
        }

        VectorFloat4 operator+(const VectorFloat4& p_vector) const {
            return VectorFloat4(
                x + p_vector.x,
                y + p_vector.y,
                z + p_vector.z,
                w + p_vector.w
            );
        }

        VectorFloat4& operator+=(const VectorFloat4& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            z += p_vector.z;
            w += p_vector.w;
            return *this;
        }

        VectorFloat4 operator-(const VectorFloat4& p_vector) const {
            return VectorFloat4(
                x - p_vector.x,
                y - p_vector.y,
                z - p_vector.z,
                w + p_vector.w
            );
        }

        VectorFloat4& operator-=(const VectorFloat4& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            z -= p_vector.z;
            w -= p_vector.w;
            return *this;
        }

        VectorFloat4 operator*(const float p_float) const {
            return VectorFloat4(
                x * p_float,
                y * p_float,
                z * p_float,
                w * p_float
            );
        }

        VectorFloat4& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            z *= p_float;
            w *= p_float;
            return *this;
        }

        VectorFloat4 operator/(const float p_float) const {
            return VectorFloat4(
                x / p_float,
                y / p_float,
                z / p_float,
                w / p_float
            );
        }

        VectorFloat4& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            z /= p_float;
            w /= p_float;
            return *this;
        }

        const float dot(const VectorFloat4& p_vector) const {
            return (
            x * p_vector.x +
            y * p_vector.y +
            z * p_vector.z +
            w * p_vector.w
            );
        }

        const float magnitude() const {
            return sqrt(
                x*x +
                y*y +
                z*z +
                w*w
            );
        }

        VectorFloat4 normalized() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return VectorFloat4(0, 0, 0, 0);
            return VectorFloat4(
                x / m,
                y / m,
                z / m,
                w / m
            );
        }

        VectorFloat4& normalize() {
            float m = magnitude();
            if (m <= FLT_EPSILON) {
                x = 0; y = 0; z = 0; w = 0;
            } else {
                x /= m;
                y /= m;
                z /= m;
                w /= m;
            }
            return *this;
        }
    };

    inline bool operator==(const VectorFloat4& p_vectorA, const VectorFloat4& p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) <= FLT_EPSILON &&
            abs(p_vectorA.y - p_vectorB.y) <= FLT_EPSILON &&
            abs(p_vectorA.z - p_vectorB.z) <= FLT_EPSILON &&
            abs(p_vectorA.w - p_vectorB.w) <= FLT_EPSILON
        );
    }

    inline bool operator!=(const VectorFloat4 &p_vectorA, const VectorFloat4 &p_vectorB) {
        return (
            abs(p_vectorA.x - p_vectorB.x) > FLT_EPSILON ||
            abs(p_vectorA.y - p_vectorB.y) > FLT_EPSILON ||
            abs(p_vectorA.z - p_vectorB.z) > FLT_EPSILON ||
            abs(p_vectorA.w - p_vectorB.w) > FLT_EPSILON
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const VectorFloat4& p_vector) {
        return p_ostream << "VectorFloat4(" << p_vector.x << ", " << p_vector.y << ", " << p_vector.z << ", " << p_vector.w << ")";
    }

} // namespace math
