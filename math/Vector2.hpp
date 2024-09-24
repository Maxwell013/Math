// Vector2.hpp
#pragma once

#include <cfloat>
#include <cmath>
#include <iostream>

namespace math {

    class Vector2 {

    public:

        float x, y;

        Vector2(const float p_x = 0, const float p_y = 0) : x(p_x), y(p_y) {}

        Vector2(const Vector2& p_vector) : x(p_vector.x), y(p_vector.y) {}

        ~Vector2() {}

        Vector2& operator=(const Vector2& p_vector) {
            x = p_vector.x;
            y = p_vector.y;
            return *this;
        }

        Vector2 operator+(const Vector2& p_vector) const {
            return Vector2(
                x + p_vector.x,
                y + p_vector.y
            );
        }

        Vector2& operator+=(const Vector2& p_vector) {
            x += p_vector.x;
            y += p_vector.y;
            return *this;
        }

        Vector2 operator-(const Vector2& p_vector) const {
            return Vector2(
                x - p_vector.x,
                y - p_vector.y
            );
        }

        Vector2& operator-=(const Vector2& p_vector) {
            x -= p_vector.x;
            y -= p_vector.y;
            return *this;
        }

        Vector2 operator*(const float p_float) const {
            return Vector2(
                x * p_float,
                y * p_float
            );
        }

        Vector2& operator*=(const float p_float) {
            x *= p_float;
            y *= p_float;
            return *this;
        }

        Vector2 operator/(const float p_float) const {
            return Vector2(
                x / p_float,
                y / p_float
            );
        }

        Vector2& operator/=(const float p_float) {
            x /= p_float;
            y /= p_float;
            return *this;
        }

        float dot(const Vector2& p_vector) const {
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

        Vector2 normalize() const {
            float m = magnitude();
            if (m <= FLT_EPSILON) return Vector2(0, 0);
            return Vector2(
                x / m,
                y / m
            );
        }

        Vector2& normalizeThis() {
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

    inline bool operator==(const Vector2& p_vector1, const Vector2& p_vector2) {
        return (
            p_vector1.x == p_vector2.x &&
            p_vector1.y == p_vector2.y
        );
    }

    inline bool operator!=(const Vector2 &p_vector1, const Vector2 &p_vector2) {
        return (
            p_vector1.x != p_vector2.x ||
            p_vector1.y != p_vector2.y
        );
    }

    inline std::ostream& operator<<(std::ostream& p_ostream, const Vector2& p_vector) {
        return p_ostream << "Vec3 : (" << p_vector.x << ", " << p_vector.y << ")!";
    }

} // namespace math
