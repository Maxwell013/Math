// math.hpp
#pragma once

// Vectors
#include "math/VectorFloat4.hpp"
#include "math/VectorFloat3.hpp"
#include "math/VectorFloat2.hpp"

// Matrices
#include "math/MatrixFloat4x4.hpp"

// Utils
#include "math/scalarFloat.hpp"
#include "math/trigonometryFloat.hpp"

namespace math {

    typedef VectorFloat4 vecf4;
    typedef VectorFloat3 vecf3;
    typedef VectorFloat2 vecf2;

    typedef vecf4 vec4;
    typedef vecf3 vec3;
    typedef vecf2 vec2;

    typedef MatrixFloat4x4 matf4x4;

    typedef matf4x4 mat4;
}
