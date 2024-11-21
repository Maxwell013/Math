// scalar.hpp
#pragma once

namespace math {

    constexpr float FLT_EPSILON = __FLT_EPSILON__;
    constexpr float FLT_MIN = __FLT_MIN__;
    constexpr float FLT_MAX = __FLT_MAX__;

    constexpr inline float min(const float p_valueA, const float p_valueB) { return (p_valueA > p_valueB) ? p_valueB : p_valueA; }

    constexpr inline float max(const float p_valueA, const float p_valueB) { return (p_valueA > p_valueB) ? p_valueA : p_valueB; }

    constexpr inline float clamp(const float p_value, const float p_min, const float p_max) {
        return min(max(p_value, p_min), p_max);
    }

    inline float abs(const float p_value) {
        unsigned int mask = 0x7FFFFFFF;
        unsigned int value = * (int *) &p_value;
        value &= mask;
        return * (float *) &value;
    }

    constexpr inline float sign(const float p_value) { return (p_value != 0) ? p_value/abs(p_value) : 0; }

    constexpr inline float lerp(const float p_value, const float p_target, const float p_fraction) {
        return (p_value * (1.0 - p_fraction)) + (p_target * p_fraction);
    }

    constexpr inline float trunc(const float p_value) {
        return (p_value > 0.0f) ? (float)((int)p_value) : (float)(-(int)(-p_value));
    }

    constexpr inline float floor(const float p_value) {
        const float value = trunc(p_value);
        return (p_value < 0 && p_value != value) ? value - 1 : value;
    }

    constexpr inline float ceil(const float p_value) {
        const float value = trunc(p_value);
        return (p_value > 0 && p_value != value) ? value + 1 : value;
    }

    constexpr inline float round(const float p_value) {
        return (p_value > 0) ? floor(p_value + 0.5f) : ceil(p_value - 0.5f);
    }

    constexpr inline float fmod(const float p_valueA, const float p_valueB) {
    float quotient = p_valueA / p_valueB;
    float integer = trunc(quotient);
    float remainder = p_valueA - integer * p_valueB;
   	return remainder + (remainder < 0.0f) * p_valueB;
    }

    constexpr inline float sqrt(const float p_value) {
    float temp;
    float root = p_value / 2;
    do {
        temp = root;
        root = (p_value/temp + temp) / 2;
    } while (root != temp);
    return root;
    }
} // namespace math
