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
      float m00, m01, m02, m03,
          m10, m11, m12, m13,
          m20, m21, m22, m23,
          m30, m31, m32, m33;
    };
    float data[16];
    vec4 vecs[4];
    //    simd_align_max float data[16];
    //    simd_float4 sse_data[4];
    //#ifdef USE_AVX2
    //    simd_float8 sse_data256[2];
    //#endif
  };
} mat4;

static inline mat4 mat4_mul(mat4 m1, mat4 m2);
static inline vec4 mat4_mul_vec4(mat4 m1, vec4 v1);
static inline mat4 mat4_translate(mat4 m1, vec3 v1);
static inline vec3 mat4_transform(mat4 m1, vec3 v1);
static inline float mat4_determinant(mat4 m1);
static inline mat4 mat4_inverse(mat4 m1);
static inline vec3 mat4_transform_direction(mat4 m1, vec3 v1);
static inline vec3 mat4_transform_inverse_direction(mat4 m1, vec3 v1);
static inline vec3 mat4_transform_inverse(mat4 m1, vec3 v1);
static inline vec3 mat4_get_axis_vector(mat4 m1, int i);
static inline mat4 mat4_orientation_and_pos(mat4 m1, quat q1, vec3 v1);
static inline void mat4_fill_gl_array(mat4 m1, float* array);
static inline mat4 mat4_look_at(vec3 eye, vec3 center, vec3 up);
static inline mat4 mat4_transpose(mat4 m1);
static inline mat4 mat4_rotate(mat4 m1, float angle, vec3 axis);

static inline mat4 mat4_mul(mat4 m1, mat4 m2) {
  //#ifdef ARCH_32_64
  //#ifdef USE_AVX2
  //  __m256 y0, y1, y2, y3, y4, y5, y6, y7, y8, y9;
  //  mat4 dest;
  //
  //  y0 = m2.sse_data256[0];
  //  y1 = m2.sse_data256[1];
  //  y2 = m1.sse_data256[0];
  //  y3 = m1.sse_data256[1];
  //  y4 = _mm256_permute2f128_ps(y2, y2, 0x03);
  //  y5 = _mm256_permute2f128_ps(y3, y3, 0x03);
  //  y6 = _mm256_permutevar_ps(y0, _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0));
  //  y7 = _mm256_permutevar_ps(y0, _mm256_set_epi32(3, 3, 3, 3, 2, 2, 2, 2));
  //  y8 = _mm256_permutevar_ps(y0, _mm256_set_epi32(0, 0, 0, 0, 1, 1, 1, 1));
  //  y9 = _mm256_permutevar_ps(y0, _mm256_set_epi32(2, 2, 2, 2, 3, 3, 3, 3));
  //
  //  dest.sse_data256[0] = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(y2, y6), _mm256_mul_ps(y3, y7)), _mm256_add_ps(_mm256_mul_ps(y4, y8), _mm256_mul_ps(y5, y9)));
  //
  //  y6 = _mm256_permutevar_ps(y1, _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0));
  //  y7 = _mm256_permutevar_ps(y1, _mm256_set_epi32(3, 3, 3, 3, 2, 2, 2, 2));
  //  y8 = _mm256_permutevar_ps(y1, _mm256_set_epi32(0, 0, 0, 0, 1, 1, 1, 1));
  //  y9 = _mm256_permutevar_ps(y1, _mm256_set_epi32(2, 2, 2, 2, 3, 3, 3, 3));
  //
  //  dest.sse_data256[1] = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(y2, y6), _mm256_mul_ps(y3, y7)), _mm256_add_ps(_mm256_mul_ps(y4, y8), _mm256_mul_ps(y5, y9)));
  //
  //  return dest;
  //#else
  //  __m128 l0, l1, l2, l3, r;
  //  mat4 dest;
  //
  //  l0 = m1.sse_data[0];
  //  l1 = m1.sse_data[1];
  //  l2 = m1.sse_data[2];
  //  l3 = m1.sse_data[3];
  //
  //  r = m2.sse_data[0];
  //  dest.sse_data[0] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  //  r = m2.sse_data[1];
  //  dest.sse_data[1] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  //  r = m2.sse_data[2];
  //  dest.sse_data[2] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  //  r = m2.sse_data[3];
  //  dest.sse_data[3] = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)), l0), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)), l1)), _mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)), l2), _mm_mul_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)), l3)));
  //
  //  return dest;
  //#endif
  //#else
  //// Neon
  //#endif
  // Reference without SIMD
  mat4 dest;
  float a00 = m1.vecs[0].data[0], a01 = m1.vecs[0].data[1], a02 = m1.vecs[0].data[2], a03 = m1.vecs[0].data[3],
        a10 = m1.vecs[1].data[0], a11 = m1.vecs[1].data[1], a12 = m1.vecs[1].data[2], a13 = m1.vecs[1].data[3],
        a20 = m1.vecs[2].data[0], a21 = m1.vecs[2].data[1], a22 = m1.vecs[2].data[2], a23 = m1.vecs[2].data[3],
        a30 = m1.vecs[3].data[0], a31 = m1.vecs[3].data[1], a32 = m1.vecs[3].data[2], a33 = m1.vecs[3].data[3],

        b00 = m2.vecs[0].data[0], b01 = m2.vecs[0].data[1], b02 = m2.vecs[0].data[2], b03 = m2.vecs[0].data[3],
        b10 = m2.vecs[1].data[0], b11 = m2.vecs[1].data[1], b12 = m2.vecs[1].data[2], b13 = m2.vecs[1].data[3],
        b20 = m2.vecs[2].data[0], b21 = m2.vecs[2].data[1], b22 = m2.vecs[2].data[2], b23 = m2.vecs[2].data[3],
        b30 = m2.vecs[3].data[0], b31 = m2.vecs[3].data[1], b32 = m2.vecs[3].data[2], b33 = m2.vecs[3].data[3];

  dest.vecs[0].data[0] = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
  dest.vecs[0].data[1] = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
  dest.vecs[0].data[2] = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
  dest.vecs[0].data[3] = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
  dest.vecs[1].data[0] = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
  dest.vecs[1].data[1] = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
  dest.vecs[1].data[2] = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
  dest.vecs[1].data[3] = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
  dest.vecs[2].data[0] = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
  dest.vecs[2].data[1] = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
  dest.vecs[2].data[2] = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
  dest.vecs[2].data[3] = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
  dest.vecs[3].data[0] = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
  dest.vecs[3].data[1] = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
  dest.vecs[3].data[2] = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
  dest.vecs[3].data[3] = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
  return dest;
  //return (mat4){.data[0] = (m2.data[0] * m1.data[0]) + (m2.data[4] * m1.data[1]) + (m2.data[8] * m1.data[2]),
  //              .data[1] = (m2.data[1] * m1.data[0]) + (m2.data[5] * m1.data[1]) + (m2.data[9] * m1.data[2]),
  //              .data[2] = (m2.data[2] * m1.data[0]) + (m2.data[6] * m1.data[1]) + (m2.data[10] * m1.data[2]),
  //              .data[3] = (m2.data[3] * m1.data[0]) + (m2.data[7] * m1.data[1]) + (m2.data[11] * m1.data[2]) + m1.data[3],
  //              .data[4] = (m2.data[0] * m1.data[4]) + (m2.data[4] * m1.data[5]) + (m2.data[8] * m1.data[6]),
  //              .data[5] = (m2.data[1] * m1.data[4]) + (m2.data[5] * m1.data[5]) + (m2.data[9] * m1.data[6]),
  //              .data[6] = (m2.data[2] * m1.data[4]) + (m2.data[6] * m1.data[5]) + (m2.data[10] * m1.data[6]),
  //              .data[7] = (m2.data[3] * m1.data[4]) + (m2.data[7] * m1.data[5]) + (m2.data[11] * m1.data[6]) + m1.data[7],
  //              .data[8] = (m2.data[0] * m1.data[8]) + (m2.data[4] * m1.data[9]) + (m2.data[8] * m1.data[10]),
  //              .data[9] = (m2.data[1] * m1.data[8]) + (m2.data[5] * m1.data[9]) + (m2.data[9] * m1.data[10]),
  //              .data[10] = (m2.data[2] * m1.data[8]) + (m2.data[6] * m1.data[9]) + (m2.data[10] * m1.data[10]),
  //              .data[11] = (m2.data[3] * m1.data[8]) + (m2.data[7] * m1.data[9]) + (m2.data[11] * m1.data[10]) + m1.data[11]};
}

static inline vec4 mat4_mul_vec4(mat4 m1, vec4 v1) {
  return (vec4){.data[0] = m1.m00 * v1.data[0] + m1.m10 * v1.data[1] + m1.m20 * v1.data[2] + m1.m30 * v1.data[3],
                .data[1] = m1.m01 * v1.data[0] + m1.m11 * v1.data[1] + m1.m21 * v1.data[2] + m1.m31 * v1.data[3],
                .data[2] = m1.m02 * v1.data[0] + m1.m12 * v1.data[1] + m1.m22 * v1.data[2] + m1.m32 * v1.data[3],
                .data[3] = m1.m03 * v1.data[0] + m1.m13 * v1.data[1] + m1.m23 * v1.data[2] + m1.m33 * v1.data[3]};
}

static inline mat4 mat4_translate(mat4 m1, vec3 v1) {
  m1.vecs[3] = vec4_add(vec4_scale(m1.vecs[0], v1.data[0]), m1.vecs[3]);
  m1.vecs[3] = vec4_add(vec4_scale(m1.vecs[1], v1.data[1]), m1.vecs[3]);
  m1.vecs[3] = vec4_add(vec4_scale(m1.vecs[2], v1.data[2]), m1.vecs[3]);
  return m1;
}

// Transform also translate? <-- NO!
static inline vec3 mat4_transform(mat4 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[1] + v1.data[2] * m1.data[2] + m1.data[3],
                .data[1] = v1.data[0] * m1.data[4] + v1.data[1] * m1.data[5] + v1.data[2] * m1.data[6] + m1.data[7],
                .data[2] = v1.data[0] * m1.data[8] + v1.data[1] * m1.data[9] + v1.data[2] * m1.data[10] + m1.data[11]};
}

static inline float mat4_determinant(mat4 m1) {
  return -m1.data[8] * m1.data[5] * m1.data[2] + m1.data[4] * m1.data[9] * m1.data[2] + m1.data[8] * m1.data[1] * m1.data[6] - m1.data[0] * m1.data[9] * m1.data[6] - m1.data[4] * m1.data[1] * m1.data[10] + m1.data[0] * m1.data[5] * m1.data[10];
}

static inline mat4 mat4_inverse(mat4 m1) {
  mat4 dest;
  float t[6];
  float det;
  float a = m1.vecs[0].data[0], b = m1.vecs[0].data[1], c = m1.vecs[0].data[2], d = m1.vecs[0].data[3],
        e = m1.vecs[1].data[0], f = m1.vecs[1].data[1], g = m1.vecs[1].data[2], h = m1.vecs[1].data[3],
        i = m1.vecs[2].data[0], j = m1.vecs[2].data[1], k = m1.vecs[2].data[2], l = m1.vecs[2].data[3],
        m = m1.vecs[3].data[0], n = m1.vecs[3].data[1], o = m1.vecs[3].data[2], p = m1.vecs[3].data[3];

  t[0] = k * p - o * l;
  t[1] = j * p - n * l;
  t[2] = j * o - n * k;
  t[3] = i * p - m * l;
  t[4] = i * o - m * k;
  t[5] = i * n - m * j;

  dest.vecs[0].data[0] = f * t[0] - g * t[1] + h * t[2];
  dest.vecs[1].data[0] = -(e * t[0] - g * t[3] + h * t[4]);
  dest.vecs[2].data[0] = e * t[1] - f * t[3] + h * t[5];
  dest.vecs[3].data[0] = -(e * t[2] - f * t[4] + g * t[5]);

  dest.vecs[0].data[1] = -(b * t[0] - c * t[1] + d * t[2]);
  dest.vecs[1].data[1] = a * t[0] - c * t[3] + d * t[4];
  dest.vecs[2].data[1] = -(a * t[1] - b * t[3] + d * t[5]);
  dest.vecs[3].data[1] = a * t[2] - b * t[4] + c * t[5];

  t[0] = g * p - o * h;
  t[1] = f * p - n * h;
  t[2] = f * o - n * g;
  t[3] = e * p - m * h;
  t[4] = e * o - m * g;
  t[5] = e * n - m * f;

  dest.vecs[0].data[2] = b * t[0] - c * t[1] + d * t[2];
  dest.vecs[1].data[2] = -(a * t[0] - c * t[3] + d * t[4]);
  dest.vecs[2].data[2] = a * t[1] - b * t[3] + d * t[5];
  dest.vecs[3].data[2] = -(a * t[2] - b * t[4] + c * t[5]);

  t[0] = g * l - k * h;
  t[1] = f * l - j * h;
  t[2] = f * k - j * g;
  t[3] = e * l - i * h;
  t[4] = e * k - i * g;
  t[5] = e * j - i * f;

  dest.vecs[0].data[3] = -(b * t[0] - c * t[1] + d * t[2]);
  dest.vecs[1].data[3] = a * t[0] - c * t[3] + d * t[4];
  dest.vecs[2].data[3] = -(a * t[1] - b * t[3] + d * t[5]);
  dest.vecs[3].data[3] = a * t[2] - b * t[4] + c * t[5];

  det = 1.0f / (a * dest.vecs[0].data[0] + b * dest.vecs[1].data[0] + c * dest.vecs[2].data[0] + d * dest.vecs[3].data[0]);

  dest.vecs[0].data[0] *= det;
  dest.vecs[0].data[1] *= det;
  dest.vecs[0].data[2] *= det;
  dest.vecs[0].data[3] *= det;
  dest.vecs[1].data[0] *= det;
  dest.vecs[1].data[1] *= det;
  dest.vecs[1].data[2] *= det;
  dest.vecs[1].data[3] *= det;
  dest.vecs[2].data[0] *= det;
  dest.vecs[2].data[1] *= det;
  dest.vecs[2].data[2] *= det;
  dest.vecs[2].data[3] *= det;
  dest.vecs[3].data[0] *= det;
  dest.vecs[3].data[1] *= det;
  dest.vecs[3].data[2] *= det;
  dest.vecs[3].data[3] *= det;
  return dest;
  //float det = mat4_determinant(m1);
  //if (det == 0)
  //  return m1;
  //det = 1.0f / det;
  //
  //return (mat4){.data[0] = (-m1.data[9] * m1.data[6] + m1.data[5] * m1.data[10]) * det,
  //              .data[1] = (m1.data[9] * m1.data[2] - m1.data[1] * m1.data[10]) * det,
  //              .data[2] = (-m1.data[5] * m1.data[2] + m1.data[1] * m1.data[6]) * det,
  //              .data[3] = (m1.data[9] * m1.data[6] * m1.data[3] - m1.data[5] * m1.data[10] * m1.data[3] - m1.data[9] * m1.data[2] * m1.data[7] + m1.data[1] * m1.data[10] * m1.data[7] + m1.data[5] * m1.data[2] * m1.data[11] - m1.data[1] * m1.data[6] * m1.data[11]) * det,
  //              .data[4] = (m1.data[8] * m1.data[6] - m1.data[4] * m1.data[10]) * det,
  //              .data[5] = (-m1.data[8] * m1.data[2] + m1.data[0] * m1.data[10]) * det,
  //              .data[6] = (+m1.data[4] * m1.data[2] - m1.data[0] * m1.data[6]) * det,
  //              .data[7] = (-m1.data[8] * m1.data[6] * m1.data[3] + m1.data[4] * m1.data[10] * m1.data[3] + m1.data[8] * m1.data[2] * m1.data[7] - m1.data[0] * m1.data[10] * m1.data[7] - m1.data[4] * m1.data[2] * m1.data[11] + m1.data[0] * m1.data[6] * m1.data[11]) * det,
  //              .data[8] = (-m1.data[8] * m1.data[5] + m1.data[4] * m1.data[9]) * det,
  //              .data[9] = (m1.data[8] * m1.data[1] - m1.data[0] * m1.data[9]) * det,
  //              .data[10] = (-m1.data[4] * m1.data[1] + m1.data[0] * m1.data[5]) * det,
  //              .data[11] = (m1.data[8] * m1.data[5] * m1.data[3] - m1.data[4] * m1.data[9] * m1.data[3] - m1.data[8] * m1.data[1] * m1.data[7] + m1.data[0] * m1.data[9] * m1.data[7] + m1.data[4] * m1.data[1] * m1.data[11] - m1.data[0] * m1.data[5] * m1.data[11]) * det};
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

static inline mat4 mat4_orientation_and_pos(mat4 m1, quat q1, vec3 v1) {
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

static inline mat4 mat4_look_at(vec3 eye, vec3 center, vec3 up) {
  mat4 dest;
  vec3 f, u, s;

  f = vec3_sub(center, eye);
  f = vec3_normalise(f);

  s = vec3_normalise(vec3_cross_product(f, up));
  u = vec3_cross_product(s, f);

  dest.vecs[0].data[0] = s.data[0];
  dest.vecs[0].data[1] = u.data[0];
  dest.vecs[0].data[2] = -f.data[0];
  dest.vecs[1].data[0] = s.data[1];
  dest.vecs[1].data[1] = u.data[1];
  dest.vecs[1].data[2] = -f.data[1];
  dest.vecs[2].data[0] = s.data[2];
  dest.vecs[2].data[1] = u.data[2];
  dest.vecs[2].data[2] = -f.data[2];
  dest.vecs[3].data[0] = -vec3_dot(s, eye);
  dest.vecs[3].data[1] = -vec3_dot(u, eye);
  dest.vecs[3].data[2] = vec3_dot(f, eye);
  dest.vecs[0].data[3] = dest.vecs[1].data[3] = dest.vecs[2].data[3] = 0.0f;
  dest.vecs[3].data[3] = 1.0f;

  return dest;
}

static inline mat4 mat4_transpose(mat4 m1) {
  //return (mat4){.m00 = m1.m00,
  //              .m10 = m1.m01,
  //              .m01 = m1.m10,
  //              .m11 = m1.m11,
  //              .m02 = m1.m20,
  //              .m12 = m1.m21,
  //              .m03 = m1.m30,
  //              .m13 = m1.m31,
  //              .m20 = m1.m02,
  //              .m30 = m1.m03,
  //              .m21 = m1.m12,
  //              .m31 = m1.m13,
  //              .m22 = m1.m22,
  //              .m32 = m1.m23,
  //              .m23 = m1.m32,
  //              .m33 = m1.m33};
  mat4 dest;
  dest.vecs[0].data[0] = m1.vecs[0].data[0];
  dest.vecs[1].data[0] = m1.vecs[0].data[1];
  dest.vecs[0].data[1] = m1.vecs[1].data[0];
  dest.vecs[1].data[1] = m1.vecs[1].data[1];
  dest.vecs[0].data[2] = m1.vecs[2].data[0];
  dest.vecs[1].data[2] = m1.vecs[2].data[1];
  dest.vecs[0].data[3] = m1.vecs[3].data[0];
  dest.vecs[1].data[3] = m1.vecs[3].data[1];
  dest.vecs[2].data[0] = m1.vecs[0].data[2];
  dest.vecs[3].data[0] = m1.vecs[0].data[3];
  dest.vecs[2].data[1] = m1.vecs[1].data[2];
  dest.vecs[3].data[1] = m1.vecs[1].data[3];
  dest.vecs[2].data[2] = m1.vecs[2].data[2];
  dest.vecs[3].data[2] = m1.vecs[2].data[3];
  dest.vecs[2].data[3] = m1.vecs[3].data[2];
  dest.vecs[3].data[3] = m1.vecs[3].data[3];
  return dest;
}

static inline mat4 mat4_rotate(mat4 m1, float angle, vec3 axis) {
  mat4 rot = MAT4_ZERO;
  vec3 axisn, v, vs;
  float c;

  c = cosf(angle);

  ////////////////////////////////////////////////////////////////
  axisn = vec3_normalise(axis);
  v = vec3_scale(axisn, 1.0f - c);
  vs = vec3_scale(axisn, sinf(angle));

  // TODO: Kinda hacky
  vec3 r1, r2, r3;
  r1 = vec3_scale(axisn, v.data[0]);
  r2 = vec3_scale(axisn, v.data[1]);
  r3 = vec3_scale(axisn, v.data[2]);

  rot.vecs[0].data[0] = r1.x;
  rot.vecs[0].data[1] = r1.y;
  rot.vecs[0].data[2] = r1.z;
  rot.vecs[1].data[0] = r2.x;
  rot.vecs[1].data[1] = r2.y;
  rot.vecs[1].data[2] = r2.z;
  rot.vecs[2].data[0] = r3.x;
  rot.vecs[2].data[1] = r3.y;
  rot.vecs[2].data[2] = r3.z;

  rot.vecs[0].data[0] += c;
  rot.vecs[1].data[0] -= vs.data[2];
  rot.vecs[2].data[0] += vs.data[1];
  rot.vecs[0].data[1] += vs.data[2];
  rot.vecs[1].data[1] += c;
  rot.vecs[2].data[1] -= vs.data[0];
  rot.vecs[0].data[2] -= vs.data[1];
  rot.vecs[1].data[2] += vs.data[0];
  rot.vecs[2].data[2] += c;

  rot.m03 = rot.m13 = rot.m23 = rot.m30 = rot.m31 = rot.m32 = 0.0f;
  rot.m33 = 1.0f;
  ////////////////////////////////////////////////////////////////
  float a00 = m1.vecs[0].data[0], a01 = m1.vecs[0].data[1], a02 = m1.vecs[0].data[2], a03 = m1.vecs[0].data[3],
        a10 = m1.vecs[1].data[0], a11 = m1.vecs[1].data[1], a12 = m1.vecs[1].data[2], a13 = m1.vecs[1].data[3],
        a20 = m1.vecs[2].data[0], a21 = m1.vecs[2].data[1], a22 = m1.vecs[2].data[2], a23 = m1.vecs[2].data[3],
        a30 = m1.vecs[3].data[0], a31 = m1.vecs[3].data[1], a32 = m1.vecs[3].data[2], a33 = m1.vecs[3].data[3],

        b00 = rot.vecs[0].data[0], b01 = rot.vecs[0].data[1], b02 = rot.vecs[0].data[2],
        b10 = rot.vecs[1].data[0], b11 = rot.vecs[1].data[1], b12 = rot.vecs[1].data[2],
        b20 = rot.vecs[2].data[0], b21 = rot.vecs[2].data[1], b22 = rot.vecs[2].data[2];

  mat4 dest = m1;
  dest.vecs[0].data[0] = a00 * b00 + a10 * b01 + a20 * b02;
  dest.vecs[0].data[1] = a01 * b00 + a11 * b01 + a21 * b02;
  dest.vecs[0].data[2] = a02 * b00 + a12 * b01 + a22 * b02;
  dest.vecs[0].data[3] = a03 * b00 + a13 * b01 + a23 * b02;
  dest.vecs[1].data[0] = a00 * b10 + a10 * b11 + a20 * b12;
  dest.vecs[1].data[1] = a01 * b10 + a11 * b11 + a21 * b12;
  dest.vecs[1].data[2] = a02 * b10 + a12 * b11 + a22 * b12;
  dest.vecs[1].data[3] = a03 * b10 + a13 * b11 + a23 * b12;
  dest.vecs[2].data[0] = a00 * b20 + a10 * b21 + a20 * b22;
  dest.vecs[2].data[1] = a01 * b20 + a11 * b21 + a21 * b22;
  dest.vecs[2].data[2] = a02 * b20 + a12 * b21 + a22 * b22;
  dest.vecs[2].data[3] = a03 * b20 + a13 * b21 + a23 * b22;
  dest.vecs[3].data[0] = a30;
  dest.vecs[3].data[1] = a31;
  dest.vecs[3].data[2] = a32;
  dest.vecs[3].data[3] = a33;
  return dest;
}

#endif  // MAT4_H
