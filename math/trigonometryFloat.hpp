// trigonometryFloat.hpp
#pragma once

#include "scalarFloat.hpp"

namespace math {

    constexpr float PI = 3.14159265358f ;
    constexpr float TWO_PI = 2.0f * PI;
    constexpr float HALF_PI = PI / 2.0f;

    constexpr inline float radians(const float p_degrees) {
        constexpr float conversion = TWO_PI / 360.0f;
        return conversion * p_degrees;
    }

    constexpr inline float degrees(const float p_radians) {
        constexpr float conversion = 360.0f / TWO_PI;
        return conversion * p_radians;
    }

    constexpr inline float sin(const float p_angle) {
    const float angle = fmod(p_angle + PI, TWO_PI) - PI;
    const float cube = angle * angle * angle;
    return angle
        -  cube                  * 0.16666666666f
        + (cube * angle * angle) * 0.00833333333f
        - (cube * cube * angle)  * 0.00019841269f
        + (cube * cube * cube)   * 0.00000275573f;
    }

    constexpr inline float cos(const float p_angle) {
    const float angle = fmod(p_angle + PI, TWO_PI) - PI;
    const float square = angle * angle;
    const float quatric = square * square;
    return 1.0f
        - square              * 0.5f
        + quatric             * 0.04166666666f
        - (quatric * square)  * 0.00138888888f
        + (quatric * quatric) * 0.00002480158f;
    }

    constexpr inline float tan(const float p_angle) {
        return sin(p_angle)/cos(p_angle);
    }
} // namespace math
