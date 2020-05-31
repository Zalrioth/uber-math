#pragma once
#ifndef VEC2_H
#define VEC2_H

#include "ubermathcommon.h"

#define VEC2_INIT_ZERO \
  { .data[0] = 0.0f, .data[1] = 0.0f }

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
  };
} vec2;

static inline vec2 vec2_add(vec2 a, vec2 b) {
  return (vec2){.data[0] = a.data[0] + b.data[0], .data[1] = a.data[1] + b.data[1]};
}

typedef struct vec2_soa {
  union {
    struct {
      float *x, *y;
    };
    struct {
      float *r, *g;
    };
    struct {
      float *u, *v;
    };
    simd_align_max float *data[2];
    struct {
      simd_float_max *simd_data_x, *simd_data_y;
    };
  };
} vec2_soa;

static inline void vec2_soa_init(vec2_soa *vec2_soa_ref) {
  vec2_soa_ref->data[0] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec2_soa_ref->data[1] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
}

static inline void vec2_soa_delete(vec2_soa *vec2_soa_ref) {
  free(vec2_soa_ref->data[0]);
  free(vec2_soa_ref->data[1]);
}

static inline void vec2_soa_resize(vec2_soa *vec2_soa_ref, size_t old_size, size_t size) {
  size_t simd_round_old = (old_size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
  size_t simd_round = (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;

  if (simd_round_old == simd_round)
    return;

  vec2_soa_ref->data[0] = realloc(vec2_soa_ref->data[0], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec2_soa_ref->data[1] = realloc(vec2_soa_ref->data[1], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
}

static inline size_t vec2_soa_iterations(size_t size) {
  return (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
}

static inline simd_float_max vec2_soa_add(simd_float_max a, simd_float_max b) {
  return simd_float_max_add(a, b);
}

#endif  // VEC2_H
