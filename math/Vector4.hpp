// Vector4.hpp
#pragma once

#include <cfloat>
#include <cmath>
#include <iostream>

namespace math {

    class Vector4 {

    public:

        float x, y, z, w;

        Vector4(const float p_x = 0, const float p_y = 0, const float p_z = 0, const float p_w = 1) : x(p_x), y(p_y), z(p_z), w(p_w) {}

        Vector4(const Vector4& p_vector) : x(p_vector.x), y(p_vector.y), z(p_vector.z), w(p_vector.w) {}

        ~Vector4() {}

        Vector4& operator=(const Vector4& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            z = p_vector.z;
            w = p_vector.w;
            return *this;
        }

        Vector4 operator+(const Vector4& p_vector) const {
            return Vector4(
                x + p_vector.x,
                y + p_vector.y,
                z + p_vector.z,
                w + p_vector.w
            );
        }

        Vector4& operator+=(const Vector4& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            z += p_vector.z;
            w += p_vector.w;
            return *this;
        }

        Vector4 operator-(const Vector4& p_vector) const {
            return Vector4(
                x - p_vector.x,
                y - p_vector.y,
                z - p_vector.z,
                w + p_vector.w
            );
        }

        Vector4& operator-=(const Vector4& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            z -= p_vector.z;
            w -= p_vector.w;
            return *this;
        }

        Vector4 operator*(const float p_float) const {
            return Vector4(
                x * p_float,
                y * p_float,
                z * p_float,
                w * p_float
            );
        }

        Vector4& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            z *= p_float;
            w *= p_float;
            return *this;
        }

        Vector4 operator/(const float p_float) const {
            return Vector4(
                x / p_float,
                y / p_float,
                z / p_float,
                w / p_float
            );
        }

        Vector4& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            z /= p_float;
            w /= p_float;
            return *this;
        }

        float dot(const Vector4& p_vector) const {
            return (
            x * p_vector.x +
            y * p_vector.y +
            z * p_vector.z +
            w * p_vector.w
            );
        }

        Vector4 cross(const Vector4& p_vector) const {
            return Vector4(
                y * p_vector.z - z * p_vector.y,
                z * p_vector.x - x * p_vector.z,
                x * p_vector.y - y * p_vector.x,
                0
            );
        }

        float magnitude() const {
            return sqrt(
                x*x +
                y*y +
                z*z +
                w*w
            );
        }

        Vector4 normalize() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return Vector4(0, 0, 0, 0);
            return Vector4(
                x / m,
                y / m,
                z / m,
                w / m
            );
        }

        Vector4& normalizeThis() {
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

    inline bool operator==(const Vector4& p_vector1, const Vector4& p_vector2) {
        return (
            p_vector1.x == p_vector2.x &&
            p_vector1.y == p_vector2.y &&
            p_vector1.z == p_vector2.z &&
            p_vector1.w == p_vector2.w
        );
    }

    inline bool operator!=(const Vector4 &p_vector1, const Vector4 &p_vector2) {
        return (
            p_vector1.x != p_vector2.x ||
            p_vector1.y != p_vector2.y ||
            p_vector1.z != p_vector2.z ||
            p_vector1.w != p_vector2.w
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const Vector4& p_vector) {
        return p_ostream << "Vec3 : (" << p_vector.x << ", " << p_vector.y << ", " << p_vector.z << ", " << p_vector.w << ")!";
    }

} // namespace math
