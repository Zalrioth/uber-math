#pragma once
#ifndef MAT4_H
#define MAT4_H

#include "quat.h"
#include "ubermathcommon.h"
#include "vec3.h"
#include "vec4.h"

#define MAT4_ZERO \
  (mat4) { .data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f, .data[4] = 0.0f, .data[5] = 0.0f, .data[6] = 0.0f, .data[7] = 0.0f, .data[8] = 0.0f, .data[9] = 0.0f, .data[10] = 0.0f, .data[11] = 0.0f, .data[12] = 0.0f, .data[13] = 0.0f, .data[14] = 0.0f, .data[15] = 0.0f }

#define MAT4_IDENTITY \
  (mat4) { .data[0] = 1.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f, .data[4] = 0.0f, .data[5] = 1.0f, .data[6] = 0.0f, .data[7] = 0.0f, .data[8] = 0.0f, .data[9] = 0.0f, .data[10] = 1.0f, .data[11] = 0.0f, .data[12] = 0.0f, .data[13] = 0.0f, .data[14] = 0.0f, .data[15] = 1.0f }

typedef struct mat4 {
  union {
    struct
    {
      float m11, m21, m31, m41,
          m12, m22, m32, m42,
          m13, m23, m33, m43,
          m14, m24, m34, m44;
    };
    vec4 vecs[4];
    simd_align_max float data[16];
    simd_float4 sse_data[4];
#ifdef USE_AVX2
    simd_float8 sse_data256[2];
#endif
  };
} mat4;

static inline mat4 mat4_mul(mat4 m1, mat4 m2);
static inline vec3 mat4_transform(mat4 m1, vec3 v1);
static inline float mat4_determinant(mat4 m1);
static inline mat4 mat4_inverse(mat4 m1);
static inline vec3 mat4_transform_direction(mat4 m1, vec3 v1);
static inline vec3 mat4_transform_inverse_direction(mat4 m1, vec3 v1);
static inline vec3 mat4_transform_inverse(mat4 m1, vec3 v1);
static inline vec3 mat4_get_axis_vector(mat4 m1, int i);
static inline mat4 mat4__orientation_and_pos(mat4 m1, quat q1, vec3 v1);
static inline void mat4_fill_gl_array(mat4 m1, float* array);

static inline mat4 mat4_mul(mat4 m1, mat4 m2) {
#ifdef ARCH_32_64
#ifdef USE_AVX2
  __m256 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9;
  mat4 dest;

  y0 = m2.sse_data256[0];
  y1 = m2.sse_data256[1];
  y2 = m1.sse_data256[0];
  y3 = m1.sse_data256[1];
  y4 = _mm256_permute2f128_ps(y2, y2, 0x03);
  y5 = _mm256_permute2f128_ps(y3, y3, 0x03);
  y6 = _mm256_permutevar_ps(y0, _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0));
  y7 = _mm256_permutevar_ps(y0, _mm256_set_epi32(3, 3, 3, 3, 2, 2, 2, 2));
  y8 = _mm256_permutevar_ps(y0, _mm256_set_epi32(0, 0, 0, 0, 1, 1, 1, 1));
  y9 = _mm256_permutevar_ps(y0, _mm256_set_epi32(2, 2, 2, 2, 3, 3, 3, 3));

  dest.sse_data256[0] = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(y2, y6), _mm256_mul_ps(y3, y7)), _mm256_add_ps(_mm256_mul_ps(y4, y8), _mm256_mul_ps(y5, y9)));

  y6 = _mm256_permutevar_ps(y1, _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0));
  y7 = _mm256_permutevar_ps(y1, _mm256_set_epi32(3, 3, 3, 3, 2, 2, 2, 2));
  y8 = _mm256_permutevar_ps(y1, _mm256_set_epi32(0, 0, 0, 0, 1, 1, 1, 1));
  y9 = _mm256_permutevar_ps(y1, _mm256_set_epi32(2, 2, 2, 2, 3, 3, 3, 3));

  dest.sse_data256[1] = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(y2, y6), _mm256_mul_ps(y3, y7)), _mm256_add_ps(_mm256_mul_ps(y4, y8), _mm256_mul_ps(y5, y9)));

  return dest;
#else
  __m128 l0, l1, l2, l3, r;
  mat4 dest;

  l0 = m1.sse_data[0];
  l1 = m1.sse_data[1];
  l2 = m1.sse_data[2];
  l3 = m1.sse_data[3];

  r = m2.sse_data[0];
  dest.sse_data[0] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  r = m2.sse_data[1];
  dest.sse_data[1] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  r = m2.sse_data[2];
  dest.sse_data[2] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  r = m2.sse_data[3];
  dest.sse_data[3] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));

  return dest;
#endif
#else
// Neon
#endif
  // Reference without SIMD
  //return (mat4){(m2.data[0] * m1.data[0]) + (m2.data[4] * m1.data[1]) + (m2.data[8] * m1.data[2]),
  //              (m2.data[1] * m1.data[0]) + (m2.data[5] * m1.data[1]) + (m2.data[9] * m1.data[2]),
  //              (m2.data[2] * m1.data[0]) + (m2.data[6] * m1.data[1]) + (m2.data[10] * m1.data[2]),
  //              (m2.data[3] * m1.data[0]) + (m2.data[7] * m1.data[1]) + (m2.data[11] * m1.data[2]) + m1.data[3],
  //              (m2.data[0] * m1.data[4]) + (m2.data[4] * m1.data[5]) + (m2.data[8] * m1.data[6]),
  //              (m2.data[1] * m1.data[4]) + (m2.data[5] * m1.data[5]) + (m2.data[9] * m1.data[6]),
  //              (m2.data[2] * m1.data[4]) + (m2.data[6] * m1.data[5]) + (m2.data[10] * m1.data[6]),
  //              (m2.data[3] * m1.data[4]) + (m2.data[7] * m1.data[5]) + (m2.data[11] * m1.data[6]) + m1.data[7],
  //              (m2.data[0] * m1.data[8]) + (m2.data[4] * m1.data[9]) + (m2.data[8] * m1.data[10]),
  //              (m2.data[1] * m1.data[8]) + (m2.data[5] * m1.data[9]) + (m2.data[9] * m1.data[10]),
  //              (m2.data[2] * m1.data[8]) + (m2.data[6] * m1.data[9]) + (m2.data[10] * m1.data[10]),
  //              (m2.data[3] * m1.data[8]) + (m2.data[7] * m1.data[9]) + (m2.data[11] * m1.data[10]) + m1.data[11]};
}

// Transform also translate?
static inline vec3 mat4_transform(mat4 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[1] + v1.data[2] * m1.data[2] + m1.data[3],
                .data[1] = v1.data[0] * m1.data[4] + v1.data[1] * m1.data[5] + v1.data[2] * m1.data[6] + m1.data[7],
                .data[2] = v1.data[0] * m1.data[8] + v1.data[1] * m1.data[9] + v1.data[2] * m1.data[10] + m1.data[11]};
}

static inline float mat4_determinant(mat4 m1) {
  return -m1.data[8] * m1.data[5] * m1.data[2] + m1.data[4] * m1.data[9] * m1.data[2] + m1.data[8] * m1.data[1] * m1.data[6] - m1.data[0] * m1.data[9] * m1.data[6] - m1.data[4] * m1.data[1] * m1.data[10] + m1.data[0] * m1.data[5] * m1.data[10];
}

static inline mat4 mat4_inverse(mat4 m1) {
  float det = mat4_determinant(m1);
  if (det == 0)
    return m1;
  det = 1.0f / det;

  return (mat4){.data[0] = (-m1.data[9] * m1.data[6] + m1.data[5] * m1.data[10]) * det,
                .data[1] = (m1.data[9] * m1.data[2] - m1.data[1] * m1.data[10]) * det,
                .data[2] = (-m1.data[5] * m1.data[2] + m1.data[1] * m1.data[6]) * det,
                .data[3] = (m1.data[9] * m1.data[6] * m1.data[3] - m1.data[5] * m1.data[10] * m1.data[3] - m1.data[9] * m1.data[2] * m1.data[7] + m1.data[1] * m1.data[10] * m1.data[7] + m1.data[5] * m1.data[2] * m1.data[11] - m1.data[1] * m1.data[6] * m1.data[11]) * det,
                .data[4] = (m1.data[8] * m1.data[6] - m1.data[4] * m1.data[10]) * det,
                .data[5] = (-m1.data[8] * m1.data[2] + m1.data[0] * m1.data[10]) * det,
                .data[6] = (+m1.data[4] * m1.data[2] - m1.data[0] * m1.data[6]) * det,
                .data[7] = (-m1.data[8] * m1.data[6] * m1.data[3] + m1.data[4] * m1.data[10] * m1.data[3] + m1.data[8] * m1.data[2] * m1.data[7] - m1.data[0] * m1.data[10] * m1.data[7] - m1.data[4] * m1.data[2] * m1.data[11] + m1.data[0] * m1.data[6] * m1.data[11]) * det,
                .data[8] = (-m1.data[8] * m1.data[5] + m1.data[4] * m1.data[9]) * det,
                .data[9] = (m1.data[8] * m1.data[1] - m1.data[0] * m1.data[9]) * det,
                .data[10] = (-m1.data[4] * m1.data[1] + m1.data[0] * m1.data[5]) * det,
                .data[11] = (m1.data[8] * m1.data[5] * m1.data[3] - m1.data[4] * m1.data[9] * m1.data[3] - m1.data[8] * m1.data[1] * m1.data[7] + m1.data[0] * m1.data[9] * m1.data[7] + m1.data[4] * m1.data[1] * m1.data[11] - m1.data[0] * m1.data[5] * m1.data[11]) * det};
}

static inline vec3 mat4_transform_direction(mat4 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[1] + v1.data[2] * m1.data[2],
                .data[1] = v1.data[0] * m1.data[4] + v1.data[1] * m1.data[5] + v1.data[2] * m1.data[6],
                .data[2] = v1.data[0] * m1.data[8] + v1.data[1] * m1.data[9] + v1.data[2] * m1.data[10]};
}

static inline vec3 mat4_transform_inverse_direction(mat4 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[4] + v1.data[2] * m1.data[8],
                .data[1] = v1.data[0] * m1.data[1] + v1.data[1] * m1.data[5] + v1.data[2] * m1.data[9],
                .data[2] = v1.data[0] * m1.data[2] + v1.data[1] * m1.data[6] + v1.data[2] * m1.data[10]};
}

static inline vec3 mat4_transform_inverse(mat4 m1, vec3 v1) {
  vec3 temp = v1;
  temp.data[0] -= m1.data[3];
  temp.data[1] -= m1.data[7];
  temp.data[2] -= m1.data[11];

  return (vec3){.data[0] = temp.data[0] * m1.data[0] + temp.data[1] * m1.data[4] + temp.data[2] * m1.data[8],
                .data[1] = temp.data[0] * m1.data[1] + temp.data[1] * m1.data[5] + temp.data[2] * m1.data[9],
                .data[2] = temp.data[0] * m1.data[2] + temp.data[1] * m1.data[6] + temp.data[2] * m1.data[10]};
}

static inline vec3 mat4_get_axis_vector(mat4 m1, int i) {
  return (vec3){.data[0] = m1.data[i], .data[1] = m1.data[i + 4], .data[2] = m1.data[i + 8]};
}

static inline mat4 mat4__orientation_and_pos(mat4 m1, quat q1, vec3 v1) {
  return (mat4){.data[0] = 1 - (2 * q1.data[2] * q1.data[2] + 2 * q1.data[3] * q1.data[3]),
                .data[1] = 2 * q1.data[1] * q1.data[2] + 2 * q1.data[3] * q1.data[0],
                .data[2] = 2 * q1.data[1] * q1.data[3] - 2 * q1.data[2] * q1.data[0],
                .data[3] = v1.data[0],
                .data[4] = 2 * q1.data[1] * q1.data[2] - 2 * q1.data[3] * q1.data[0],
                .data[5] = 1 - (2 * q1.data[1] * q1.data[1] + 2 * q1.data[3] * q1.data[3]),
                .data[6] = 2 * q1.data[2] * q1.data[3] + 2 * q1.data[1] * q1.data[0],
                .data[7] = v1.data[1],
                .data[8] = 2 * q1.data[1] * q1.data[3] + 2 * q1.data[2] * q1.data[0],
                .data[9] = 2 * q1.data[2] * q1.data[3] - 2 * q1.data[1] * q1.data[0],
                .data[10] = 1 - (2 * q1.data[1] * q1.data[1] + 2 * q1.data[2] * q1.data[2]),
                .data[11] = v1.data[2]};
}

static inline void mat4_fill_gl_array(mat4 m1, float* array) {
  array[0] = m1.data[0];
  array[1] = m1.data[4];
  array[2] = m1.data[8];
  array[3] = 0.0f;

  array[4] = m1.data[1];
  array[5] = m1.data[5];
  array[6] = m1.data[9];
  array[7] = 0.0f;

  array[8] = m1.data[2];
  array[9] = m1.data[6];
  array[10] = m1.data[10];
  array[11] = 0.0f;

  array[12] = m1.data[3];
  array[13] = m1.data[7];
  array[14] = m1.data[11];
  array[15] = 1.0f;
}

#endif  // MAT4_H
