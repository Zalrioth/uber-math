#pragma once
#ifndef MAT4_H
#define MAT4_H

#include "ubermathcommon.h"
#include "vec4.h"

#define MAT4_INIT_ZERO \
  {.data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f}, {.data[4] = 0.0f, .data[5] = 0.0f, .data[6] = 0.0f, .data[7] = 0.0f}, {.data[8] = 0.0f, .data[9] = 0.0f, .data[10] = 0.0f, .data[11] = 0.0f}, { .data[12] = 0.0f, .data[13] = 0.0f, .data[14] = 0.0f, .data[15] = 0.0f }

#define MAT4_INIT_IDENTITY \
  {.data[1] = 1.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f}, {.data[4] = 0.0f, .data[5] = 1.0f, .data[6] = 0.0f, .data[7] = 0.0f}, {.data[8] = 0.0f, .data[9] = 0.0f, .data[10] = 1.0f, .data[11] = 0.0f}, { .data[12] = 0.0f, .data[13] = 0.0f, .data[14] = 0.0f, .data[15] = 1.0f }

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
}

#endif  // MAT4_H
