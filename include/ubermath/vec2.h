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
    simd_align_max float data[2];
  };
} vec2;

static inline vec2 vec2_add(vec2 a, vec2 b) {
  return (vec2){.data[0] = a.data[0] + b.data[0], .data[1] = a.data[1] + b.data[1]};
}

static inline vec2 vec2_sub(vec2 a, vec2 b) {
  return (vec2){.data[0] = a.data[0] - b.data[0], .data[1] = a.data[1] - b.data[1]};
}

static inline vec2 vec2_mul(vec2 a, vec2 b) {
  return (vec2){.data[0] = a.data[0] * b.data[0], .data[1] = a.data[1] * b.data[1]};
}

static inline vec2 vec2_div(vec2 a, vec2 b) {
  return (vec2){.data[0] = a.data[0] / b.data[0], .data[1] = a.data[1] / b.data[1]};
}

static inline vec2 vec2_divs(vec2 a, float s) {
  return (vec2){.data[0] = a.data[0] / s, .data[1] = a.data[1] / s};
}

static inline vec2 vec2_scale(vec2 a, float s) {
  return (vec2){.data[0] = a.data[0] * s, .data[1] = a.data[1] * s};
}

static inline float vec2_dot(vec2 a, vec2 b) {
  return a.data[0] * b.data[0] + a.data[1] * b.data[1];
}

typedef struct vec2_simd {
  union {
    struct {
      float x[SIMD_MAX_LENGTH], y[SIMD_MAX_LENGTH];
    };
    struct {
      float r[SIMD_MAX_LENGTH], g[SIMD_MAX_LENGTH];
    };
    struct {
      float u[SIMD_MAX_LENGTH], v[SIMD_MAX_LENGTH];
    };
    simd_align_max float data[2][SIMD_MAX_LENGTH];
    struct {
      simd_float_max simd_data[2];
    };
  };
} vec2_simd;

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
      simd_float_max *simd_data[2];
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

static inline void vec2_soa_add(vec2_soa *a, size_t iter_num, vec2_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_add(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_add(a->simd_data[1][iter_num], b.simd_data[1]);
}

static inline void vec2_soa_sub(vec2_soa *a, size_t iter_num, vec2_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_sub(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_sub(a->simd_data[1][iter_num], b.simd_data[1]);
}

static inline void vec2_soa_mul(vec2_soa *a, size_t iter_num, vec2_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1]);
}

static inline void vec2_soa_div(vec2_soa *a, size_t iter_num, vec2_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], b.simd_data[1]);
}

static inline void vec2_soa_divs(vec2_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], simd_float_max_set1(s));
}

static inline void vec2_soa_scale(vec2_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], simd_float_max_set1(s));
}

static inline simd_float_max vec2_soa_dot(vec2_soa *a, size_t iter_num, vec2_simd b) {
  return simd_float_max_add(simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]), simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1]));
}

#endif  // VEC2_H
