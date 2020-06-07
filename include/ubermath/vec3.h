#pragma once
#ifndef VEC3_H
#define VEC3_H

#include "ubermathcommon.h"

#define VEC3_INIT_ZERO \
  { .data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f }

typedef struct vec3 {
  union {
    struct {
      float x, y, z;
    };
    struct {
      float r, g, b;
    };
    struct {
      float u, v, s;
    };
    float data[3];
  };
} vec3;

static inline vec3 vec3_add(vec3 a, vec3 b) {
  return (vec3){.data[0] = a.data[0] + b.data[0], .data[1] = a.data[1] + b.data[1], .data[2] = a.data[2] + b.data[2]};
}

static inline vec3 vec3_sub(vec3 a, vec3 b) {
  return (vec3){.data[0] = a.data[0] - b.data[0], .data[1] = a.data[1] - b.data[1], .data[2] = a.data[2] - b.data[2]};
}

static inline vec3 vec3_mul(vec3 a, vec3 b) {
  return (vec3){.data[0] = a.data[0] * b.data[0], .data[1] = a.data[1] * b.data[1], .data[2] = a.data[2] * b.data[2]};
}

static inline vec3 vec3_div(vec3 a, vec3 b) {
  return (vec3){.data[0] = a.data[0] / b.data[0], .data[1] = a.data[1] / b.data[1], .data[2] = a.data[2] / b.data[2]};
}

static inline vec3 vec3_divs(vec3 a, float s) {
  return (vec3){.data[0] = a.data[0] / s, .data[1] = a.data[1] / s, .data[2] = a.data[2] / s};
}

static inline vec3 vec3_scale(vec3 a, float s) {
  return (vec3){.data[0] = a.data[0] * s, .data[1] = a.data[1] * s, .data[2] = a.data[2] * s};
}

static inline float vec3_dot(vec3 a, vec3 b) {
  return a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2];
}

typedef struct vec3_simd {
  union {
    struct {
      float x[SIMD_MAX_LENGTH], y[SIMD_MAX_LENGTH], z[SIMD_MAX_LENGTH];
    };
    struct {
      float r[SIMD_MAX_LENGTH], g[SIMD_MAX_LENGTH], b[SIMD_MAX_LENGTH];
    };
    struct {
      float u[SIMD_MAX_LENGTH], v[SIMD_MAX_LENGTH], s[SIMD_MAX_LENGTH];
    };
    simd_align_max float data[3][SIMD_MAX_LENGTH];
    struct {
      simd_float_max simd_data[3];
    };
  };
} vec3_simd;

typedef struct vec3_soa {
  union {
    struct {
      float *x, *y, *z;
    };
    struct {
      float *r, *g, *b;
    };
    struct {
      float *u, *v, *s;
    };
    simd_align_max float *data[3];
    struct {
      simd_float_max *simd_data[3];
    };
  };
} vec3_soa;

static inline void vec3_soa_init(vec3_soa *vec3_soa_ref) {
  vec3_soa_ref->data[0] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec3_soa_ref->data[1] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
  vec3_soa_ref->data[2] = malloc(sizeof(float) * SIMD_MAX_LENGTH);
}

static inline void vec3_soa_delete(vec3_soa *vec3_soa_ref) {
  free(vec3_soa_ref->data[0]);
  free(vec3_soa_ref->data[1]);
  free(vec3_soa_ref->data[2]);
}

static inline void vec3_soa_resize(vec3_soa *vec3_soa_ref, size_t old_size, size_t size) {
  size_t simd_round_old = (old_size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
  size_t simd_round = (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;

  if (simd_round_old == simd_round)
    return;

  vec3_soa_ref->data[0] = realloc(vec3_soa_ref->data[0], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec3_soa_ref->data[1] = realloc(vec3_soa_ref->data[1], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
  vec3_soa_ref->data[2] = realloc(vec3_soa_ref->data[2], sizeof(float) * (SIMD_MAX_LENGTH * simd_round));
}

static inline size_t vec3_soa_iterations(size_t size) {
  return (size + SIMD_MAX_LENGTH - 1) / SIMD_MAX_LENGTH;
}

static inline void vec3_soa_add(vec3_soa *a, size_t iter_num, vec3_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_add(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_add(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_add(a->simd_data[2][iter_num], b.simd_data[2]);
}

static inline void vec3_soa_sub(vec3_soa *a, size_t iter_num, vec3_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_sub(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_sub(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_sub(a->simd_data[2][iter_num], b.simd_data[2]);
}

static inline void vec3_soa_mul(vec3_soa *a, size_t iter_num, vec3_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_mul(a->simd_data[2][iter_num], b.simd_data[2]);
}

static inline void vec3_soa_div(vec3_soa *a, size_t iter_num, vec3_simd b) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], b.simd_data[0]);
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], b.simd_data[1]);
  a->simd_data[2][iter_num] = simd_float_max_div(a->simd_data[2][iter_num], b.simd_data[2]);
}

static inline void vec3_soa_divs(vec3_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_div(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_div(a->simd_data[1][iter_num], simd_float_max_set1(s));
  a->simd_data[2][iter_num] = simd_float_max_div(a->simd_data[2][iter_num], simd_float_max_set1(s));
}

static inline void vec3_soa_scale(vec3_soa *a, size_t iter_num, float s) {
  a->simd_data[0][iter_num] = simd_float_max_mul(a->simd_data[0][iter_num], simd_float_max_set1(s));
  a->simd_data[1][iter_num] = simd_float_max_mul(a->simd_data[1][iter_num], simd_float_max_set1(s));
  a->simd_data[2][iter_num] = simd_float_max_mul(a->simd_data[2][iter_num], simd_float_max_set1(s));
}

static inline simd_float_max vec3_soa_dot(vec3_soa *a, size_t iter_num, vec3_simd b) {
  return simd_float_max_add(simd_float_max_add(simd_float_max_mul(a->simd_data[0][iter_num], b.simd_data[0]), simd_float_max_mul(a->simd_data[1][iter_num], b.simd_data[1])), simd_float_max_mul(a->simd_data[2][iter_num], b.simd_data[2]));
}

#endif  // VEC3_H
