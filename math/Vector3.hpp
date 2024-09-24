// Vector3.hpp
#pragma once

#include <cfloat>
#include <cmath>
#include <iostream>

namespace math {

    class Vector3 {

    public:

        float x, y, z;

        Vector3(const float p_x = 0, const float p_y = 0, const float p_z = 0) : x(p_x), y(p_y), z(p_z) {}

        Vector3(const Vector3& p_vector) : x(p_vector.x), y(p_vector.y), z(p_vector.z) {}

        ~Vector3() {}

        Vector3& operator=(const Vector3& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            z = p_vector.z;
            return *this;
        }

        Vector3 operator+(const Vector3& p_vector) const {
            return Vector3(
                x + p_vector.x,
                y + p_vector.y,
                z + p_vector.z
            );
        }

        Vector3& operator+=(const Vector3& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            z += p_vector.z;
            return *this;
        }

        Vector3 operator-(const Vector3& p_vector) const {
            return Vector3(
                x - p_vector.x,
                y - p_vector.y,
                z - p_vector.z
            );
        }

        Vector3& operator-=(const Vector3& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            z -= p_vector.z;
            return *this;
        }

        Vector3 operator*(const float p_float) const {
            return Vector3(
                x * p_float,
                y * p_float,
                z * p_float
            );
        }

        Vector3& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            z *= p_float;
            return *this;
        }

        Vector3 operator/(const float p_float) const {
            return Vector3(
                x / p_float,
                y / p_float,
                z / p_float
            );
        }

        Vector3& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            z /= p_float;
            return *this;
        }

        float dot(const Vector3& p_vector) const {
            return (
            x * p_vector.x +
            y * p_vector.y +
            z * p_vector.z
            );
        }

        Vector3 cross(const Vector3& p_vector) const {
            return Vector3(
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

        Vector3 normalize() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return Vector3(0, 0, 0);
            return Vector3(
                x / m,
                y / m,
                z / m
            );
        }

        Vector3& normalizeThis() {
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

    inline bool operator==(const Vector3& p_vector1, const Vector3& p_vector2) {
        return (
            p_vector1.x == p_vector2.x &&
            p_vector1.y == p_vector2.y &&
            p_vector1.z == p_vector2.z
        );
    }

    inline bool operator!=(const Vector3 &p_vector1, const Vector3 &p_vector2) {
        return (
            p_vector1.x != p_vector2.x ||
            p_vector1.y != p_vector2.y ||
            p_vector1.z != p_vector2.z
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const Vector3& p_vector) {
        return p_ostream << "Vec3 : (" << p_vector.x << ", " << p_vector.y << ", " << p_vector.z << ")!";
    }

} // namespace math
