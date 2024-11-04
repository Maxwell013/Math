#pragma once
// simd.hpp

#if defined (__ARM_NEON) || defined (__ARM_NEON__)
    #include <arm_neon.h>
    #define SIMD_ARCH_ARM_NEON
#elif defined (__AVX__)
    #include <immintrin.h>
    #define SIMD_ARCH_X86_AVX
#elif defined (__SSE__)
    #include <xmmintrin.h>
    #define SIMD_ARCH_X86_SSE
#else
    #define SIMD_ARCH_NONE
#endif

#ifdef SIMD_ARCH_ARM_NEON
    #define SIMD_F32X4_T float32x4_t

    #define SIMD_LOAD_F32(array) vld1q_f32(array)
    #define SIMD_DUP_F32(array) vdupq_n_f32(array)
    #define SIMD_MOVE_F32(value) vmovq_n_f32(value)
    #define SIMD_STORE_F32(array, register) vst1q_f32(array, register)

    #define SIMD_ADD_F32(register1, register2) vaddq_f32(register1, register2)
    #define SIMD_REC_F32(value) vrecpeq_f32(value)
    #define SIMD_MUL_F32(register1, register2) vmulq_f32(register1, register2)

    #define SIMD_FMA_F32(register1, register2, register3, lane) vfmaq_laneq_f32(register1, register2, register3, lane)
    #define SIMD_MLA_S_F32(register1, register2, register3, lane) vmlaq_laneq_f32(register1, register2, register3)
#endif
