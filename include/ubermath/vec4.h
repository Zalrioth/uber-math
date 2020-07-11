#pragma once
#ifndef VEC4_H
#define VEC4_H

#include "ubermathcommon.h"

#define VEC4_ZERO \
  (vec4) { .data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f }

typedef struct vec4 {
  union {
    struct {
      float x, y, z, w;
    };
    struct {
      float r, g, b, a;
    };
    struct {
      float u, v, s, t;
    };
    simd_align_max float data[4];
  };
} vec4;

static inline vec4 vec4_add(vec4 a, vec4 b) {
  return (vec4){.data[0] = a.data[0] + b.data[0], .data[1] = a.data[1] + b.data[1], .data[2] = a.data[2] + b.data[2], .data[3] = a.data[3] + b.data[3]};
}

static inline vec4 vec4_sub(vec4 a, vec4 b) {
  return (vec4){.data[0] = a.data[0] - b.data[0], .data[1] = a.data[1] - b.data[1], .data[2] = a.data[2] - b.data[2], .data[3] = a.data[3] - b.data[3]};
}

static inline vec4 vec4_mul(vec4 a, vec4 b) {
  return (vec4){.data[0] = a.data[0] * b.data[0], .data[1] = a.data[1] * b.data[1], .data[2] = a.data[2] * b.data[2], .data[3] = a.data[3] * b.data[3]};
}

static inline vec4 vec4_div(vec4 a, vec4 b) {
  return (vec4){.data[0] = a.data[0] / b.data[0], .data[1] = a.data[1] / b.data[1], .data[2] = a.data[2] / b.data[2], .data[3] = a.data[3] / b.data[3]};
}

static inline vec4 vec4_divs(vec4 a, float s) {
  return (vec4){.data[0] = a.data[0] / s, .data[1] = a.data[1] / s, .data[2] = a.data[2] / s, .data[3] = a.data[3] / s};
}

static inline vec4 vec4_scale(vec4 a, float s) {
  return (vec4){.data[0] = a.data[0] * s, .data[1] = a.data[1] * s, .data[2] = a.data[2] * s, .data[3] = a.data[3] * s};
}

static inline float vec4_dot(vec4 a, vec4 b) {
  return a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2] + a.data[3] * b.data[3];
}

typedef struct vec4_simd {
  union {
    struct {
      float x[SIMD_MAX_LENGTH], y[SIMD_MAX_LENGTH], z[SIMD_MAX_LENGTH], w[SIMD_MAX_LENGTH];
    };
    struct {
      float r[SIMD_MAX_LENGTH], g[SIMD_MAX_LENGTH], b[SIMD_MAX_LENGTH], a[SIMD_MAX_LENGTH];
    };
    struct {
      float u[SIMD_MAX_LENGTH], v[SIMD_MAX_LENGTH], s[SIMD_MAX_LENGTH], t[SIMD_MAX_LENGTH];
    };
    simd_align_max float data[4][SIMD_MAX_LENGTH];
    struct {
      simd_float_max simd_data[4];
    };
  };
} vec4_simd;

typedef struct vec4_soa {
  union {
    struct {
      float *x, *y, *z, *w;
    };
    struct {
      float *r, *g, *b, *a;
    };
    struct {
      float *u, *v, *s, *t;
    };
    simd_align_max float *data[4];
    struct {
      simd_float_max *simd_data[4];
    };
  };
} vec4_soa;

static inline void vec4_soa_init(vec4_soa *vec4_soa_ref) {
  vec4_soa_ref->data[0] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec4_soa_ref->data[1] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec4_soa_ref->data[2] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec4_soa_ref->data[3] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
}

static inline void vec4_soa_delete(vec4_soa *vec4_soa_ref) {
  free(vec4_soa_ref->data[0]);
  free(vec4_soa_ref->data[1]);
  free(vec4_soa_ref->data[2]);
  free(vec4_soa_ref->data[3]);
}

static inline void vec4_soa_resize(vec4_soa *vec4_soa_ref, size_t old_size, size_t size) {
  size_t simd_round_old = (old_size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
  size_t simd_round = (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;

  if (simd_round_old == simd_round)
    return;

  vec4_soa_ref->data[0] = realloc(vec4_soa_ref->data[0], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec4_soa_ref->data[1] = realloc(vec4_soa_ref->data[1], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec4_soa_ref->data[2] = realloc(vec4_soa_ref->data[2], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec4_soa_ref->data[3] = realloc(vec4_soa_ref->data[2], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
}

static inline size_t vec4_soa_iterations(size_t size) {
  return (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
}

static inline void vec4_soa_add(vec4_soa *a, size_t iter_num, vec4_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_add(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_add(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_add(a->simd_data[2][iter_num], b.simd_data[2]);
  a->simd_data[3][iter_num] = simd_float_max_add(a->simd_data[3][iter_num], b.simd_data[3]);
}

static inline void vec4_soa_sub(vec4_soa *a, size_t iter_num, vec4_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_sub(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_sub(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_sub(a->simd_data[2][iter_num], b.simd_data[2]);
  a->simd_data[3][iter_num] = simd_float_max_sub(a->simd_data[3][iter_num], b.simd_data[3]);
}

static inline void vec4_soa_mul(vec4_soa *a, size_t iter_num, vec4_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_mul(a->simd_data[2][iter_num], b.simd_data[2]);
  a->simd_data[3][iter_num] = simd_float_max_mul(a->simd_data[3][iter_num], b.simd_data[3]);
}

static inline void vec4_soa_div(vec4_soa *a, size_t iter_num, vec4_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_div(a->simd_data[2][iter_num], b.simd_data[2]);
  a->simd_data[3][iter_num] = simd_float_max_div(a->simd_data[3][iter_num], b.simd_data[3]);
}

static inline void vec4_soa_divs(vec4_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], simd_float_max_set1(s));
  a->simd_data[2][iter_num] = simd_float_max_div(a->simd_data[2][iter_num], simd_float_max_set1(s));
  a->simd_data[3][iter_num] = simd_float_max_div(a->simd_data[3][iter_num], simd_float_max_set1(s));
}

static inline void vec4_soa_scale(vec4_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], simd_float_max_set1(s));
  a->simd_data[2][iter_num] = simd_float_max_mul(a->simd_data[2][iter_num], simd_float_max_set1(s));
  a->simd_data[3][iter_num] = simd_float_max_mul(a->simd_data[3][iter_num], simd_float_max_set1(s));
}

static inline simd_float_max vec4_soa_dot(vec4_soa *a, size_t iter_num, vec4_simd b) {
  return simd_float_max_add(simd_float_max_add(simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]), simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1])), simd_float_max_add(simd_float_max_mul(a->simd_data[2][iter_num], b.simd_data[2]), simd_float_max_mul(a->simd_data[3][iter_num], b.simd_data[3])));
}

#endif  // VEC4_H
