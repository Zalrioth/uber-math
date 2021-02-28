#pragma once
#ifndef UBER_MATH_COMMON_H
#define UBER_MATH_COMMON_H

#include <math.h>
#include <omp.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MATH_SIMD
#if defined(__arm__) || defined(__aarch64__)
#define ARCH_ARM
#include <arm_neon.h>
#else
#define ARCH_32_64
#include <smmintrin.h>

// TODO: This should be modified by user to toggle
#define USE_AVX2

#define simd_align16 alignas(16)
#define simd_float4 __m128
#define simd_float4_set1 _mm_set1_ps

#ifdef USE_AVX2
#include <immintrin.h>
#define SIMD_MAX_LENGTH 8

#define simd_align32 alignas(32)
#define simd_float8 __m256
#define simd_float8_set1 _mm256_set1_ps

#define simd_align_max simd_align32
#define simd_float_max simd_float8
#define simd_float_max_set1 simd_float8_set1
#define simd_float_max_add _mm256_add_ps
#define simd_float_max_sub _mm256_sub_ps
#define simd_float_max_mul _mm256_mul_ps
#define simd_float_max_div _mm256_div_ps
#else
#define SIMD_MAX_LENGTH 4

#define simd_align_max simd_align16
#define simd_float_max simd_float4
#define simd_float_max_set1 simd_float4_set1
#define simd_float_max_add _mm_add_ps
#define simd_float_max_sub _mm_sub_ps
#define simd_float_max_mul _mm_mul_ps
#define simd_float_max_div _mm_div_ps
#endif
#endif
#endif

#endif  // UBER_MATH_COMMON_H
