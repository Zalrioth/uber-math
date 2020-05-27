#pragma once
#ifndef VEC2_H
#define VEC2_H

#include "ubermathcommon.h"

#define VEC2_INIT_ZERO \
  { .x = 0.0f, .y = 0.0f }

// TODO: Probably an automatic way to do this with macros
#ifdef USE_AVX2
#define VEC2_SOA_INIT_ZERO \
  { .x[0] = 0.0f, .y[0] = 0.0f, x[1] = 0.0f, .y[1] = 0.0f, x[2] = 0.0f, .y[2] = 0.0f, x[3] = 0.0f, .y[3] = 0.0f, x[4] = 0.0f, .y[4] = 0.0f, x[5] = 0.0f, .y[5] = 0.0f, x[6] = 0.0f, .y[6] = 0.0f, x[7] = 0.0f, .y[7] = 0.0f }
#else
#define VEC2_SOA_INIT_ZERO \
  { .x[0] = 0.0f, .y[0] = 0.0f, x[1] = 0.0f, .y[1] = 0.0f, x[2] = 0.0f, .y[2] = 0.0f, x[3] = 0.0f, .y[3] = 0.0f, x[4] = 0.0f, .y[4] = 0.0f }
#endif

typedef struct vec2 {
  union {
    struct {
      float x, y;
    };
    struct {
      float r, g;
    };
    struct {
      float u, v;
    };
    float data[2];
    simd_float4 simd_data;
  };
} vec2;

static inline vec2 vec2_add(vec2 a, vec2 b) {
  return (vec2){.simd_data = _mm_add_ps(a.simd_data, b.simd_data)};
}

typedef struct vec2_soa {
  union {
    struct {
      float x[SIMD_LENGTH], y[SIMD_LENGTH];
    };
    struct {
      float r[SIMD_LENGTH], g[SIMD_LENGTH];
    };
    struct {
      float u[SIMD_LENGTH], v[SIMD_LENGTH];
    };
    float data[2][SIMD_LENGTH];
    struct {
      simd_float_max simd_data_x;
      simd_float_max simd_data_y;
    };
  };
} vec2_soa;

static inline vec2_soa vec2_soa_add(vec2_soa a, vec2_soa b) {
  return (vec2_soa){.simd_data_x = _mm_add_ps(a.simd_data_x, b.simd_data_x), .simd_data_y = _mm_add_ps(a.simd_data_y, b.simd_data_y)};
}

#endif  // VEC2_H
