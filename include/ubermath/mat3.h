#pragma once
#ifndef MAT3_H
#define MAT3_H

#include "quat.h"
#include "ubermathcommon.h"
#include "vec3.h"

#define MAT3_ZERO \
  (mat3) { .data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f, .data[4] = 0.0f, .data[5] = 0.0f, .data[6] = 0.0f, .data[7] = 0.0f, .data[8] = 0.0f }

#define MAT3_IDENTITY \
  (mat3) { .data[0] = 1.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f, .data[4] = 1.0f, .data[5] = 0.0f, .data[6] = 0.0f, .data[7] = 0.0f, .data[8] = 1.0f }

typedef struct mat3 {
  union {
    struct
    {
      float m00, m01, m02,
          m10, m11, m12,
          m20, m21, m22;
    };
    float data[9];
    vec3 vecs[3];
    //simd_align_max float data[9];
  };
} mat3;

static inline mat3 mat_inertia_tensor_coeffs(float ix, float iy, float iz, float ixy, float ixz, float iyz);
static inline mat3 mat3_block_intertia_tensor(mat3 m1, vec3 v1, float mass);
static inline mat3 mat3_skew_symmetric(mat3 m1, vec3 v1);
static inline vec3 mat3_transform(mat3 m1, vec3 v1);
static inline vec3 mat3_transform_transpose(mat3 m1, vec3 v1);
static inline vec3 mat3_get_row_vector(mat3 m1, int i);
static inline vec3 mat3_get_axis_vector(mat3 m1, int i);
static inline mat3 mat3_inverse(mat3 m1);
static inline mat3 mat3_transpose(mat3 m1);
static inline mat3 mat3_mul_mat3(mat3 m1, mat3 m2);
static inline mat3 mat3_mul_scalar(mat3 m1, float scalar);
static inline mat3 mat3_add_mat3(mat3 m1, mat3 m2);
static inline mat3 mat3_orientation(mat3 q1);
static inline mat3 mat3_liner_interpolate(mat3 m1, mat3 m2, float prop);

static inline mat3 mat3_block_intertia_tensor(mat3 m1, vec3 v1, float mass) {
  vec3 squares = vec3_component_product(v1, v1);
  vec3 tensor_val = (vec3){.x = 0.3f * mass * (squares.data[1] + squares.data[2]),
                           .y = 0.3f * mass * (squares.data[0] + squares.data[2]),
                           .z = 0.3f * mass * (squares.data[0] + squares.data[1])};

  return (mat3){.data[0] = tensor_val.x,
                .data[1] = 0.0f,
                .data[2] = 0.0f,
                .data[3] = 0.0f,
                .data[4] = tensor_val.y,
                .data[5] = 0.0f,
                .data[6] = 0.0f,
                .data[7] = 0.0f,
                .data[8] = tensor_val.z};
}

static inline mat3 mat3_skew_symmetric(mat3 m1, vec3 v1) {
  return (mat3){.data[0] = 0,
                .data[1] = -v1.data[2],
                .data[2] = v1.data[1],
                .data[3] = v1.data[2],
                .data[4] = 0,
                .data[5] = -v1.data[0],
                .data[6] = -v1.data[1],
                .data[7] = v1.data[0],
                .data[8] = 0};
}

// Mul vector
static inline vec3 mat3_transform(mat3 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[1] + v1.data[2] * m1.data[2],
                .data[1] = v1.data[0] * m1.data[3] + v1.data[1] * m1.data[4] + v1.data[2] * m1.data[5],
                .data[2] = v1.data[0] * m1.data[6] + v1.data[1] * m1.data[7] + v1.data[2] * m1.data[8]};
}

static inline vec3 mat3_transform_transpose(mat3 m1, vec3 v1) {
  return (vec3){.data[0] = v1.data[0] * m1.data[0] + v1.data[1] * m1.data[3] + v1.data[2] * m1.data[6],
                .data[1] = v1.data[0] * m1.data[1] + v1.data[1] * m1.data[4] + v1.data[2] * m1.data[7],
                .data[2] = v1.data[0] * m1.data[2] + v1.data[1] * m1.data[5] + v1.data[2] * m1.data[8]};
}

static inline vec3 mat3_get_row_vector(mat3 m1, int i) {
  return (vec3){.data[0] = m1.data[i * 3], .data[1] = m1.data[i * 3 + 1], .data[2] = m1.data[i * 3 + 2]};
}

static inline vec3 mat3_get_axis_vector(mat3 m1, int i) {
  return (vec3){.data[0] = m1.data[i], .data[1] = m1.data[i + 3], .data[2] = m1.data[i + 6]};
}

static inline mat3 mat3_inverse(mat3 m1) {
  float t4 = m1.data[0] * m1.data[4];
  float t6 = m1.data[0] * m1.data[5];
  float t8 = m1.data[1] * m1.data[3];
  float t10 = m1.data[2] * m1.data[3];
  float t12 = m1.data[1] * m1.data[6];
  float t14 = m1.data[2] * m1.data[6];
  float t16 = (t4 * m1.data[8] - t6 * m1.data[7] - t8 * m1.data[8] + t10 * m1.data[7] + t12 * m1.data[5] - t14 * m1.data[4]);
  if (t16 == 0.0f)
    return m1;

  float t17 = 1.0f / t16;
  return (mat3){.data[0] = (m1.data[4] * m1.data[8] - m1.data[5] * m1.data[7]) * t17,
                .data[1] = -(m1.data[1] * m1.data[8] - m1.data[2] * m1.data[7]) * t17,
                .data[2] = (m1.data[1] * m1.data[5] - m1.data[2] * m1.data[4]) * t17,
                .data[3] = -(m1.data[3] * m1.data[8] - m1.data[5] * m1.data[6]) * t17,
                .data[4] = (m1.data[0] * m1.data[8] - t14) * t17,
                .data[5] = -(t6 - t10) * t17,
                .data[6] = (m1.data[3] * m1.data[7] - m1.data[4] * m1.data[6]) * t17,
                .data[7] = -(m1.data[0] * m1.data[7] - t12) * t17,
                .data[8] = (t4 - t8) * t17};
}

static inline mat3 mat3_transpose(mat3 m1) {
  return (mat3){.data[0] = m1.data[0],
                .data[1] = m1.data[3],
                .data[2] = m1.data[6],
                .data[3] = m1.data[1],
                .data[4] = m1.data[4],
                .data[5] = m1.data[7],
                .data[6] = m1.data[2],
                .data[7] = m1.data[5],
                .data[8] = m1.data[8]};
}

static inline mat3 mat3_mul_mat3(mat3 m1, mat3 m2) {
  return (mat3){.data[0] = m1.data[0] * m2.data[0] + m1.data[1] * m2.data[3] + m1.data[2] * m2.data[6],
                .data[1] = m1.data[0] * m2.data[1] + m1.data[1] * m2.data[4] + m1.data[2] * m2.data[7],
                .data[2] = m1.data[0] * m2.data[2] + m1.data[1] * m2.data[5] + m1.data[2] * m2.data[8],
                .data[3] = m1.data[3] * m2.data[0] + m1.data[4] * m2.data[3] + m1.data[5] * m2.data[6],
                .data[4] = m1.data[3] * m2.data[1] + m1.data[4] * m2.data[4] + m1.data[5] * m2.data[7],
                .data[5] = m1.data[3] * m2.data[2] + m1.data[4] * m2.data[5] + m1.data[5] * m2.data[8],
                .data[6] = m1.data[6] * m2.data[0] + m1.data[7] * m2.data[3] + m1.data[8] * m2.data[6],
                .data[7] = m1.data[6] * m2.data[1] + m1.data[7] * m2.data[4] + m1.data[8] * m2.data[7],
                .data[8] = m1.data[6] * m2.data[2] + m1.data[7] * m2.data[5] + m1.data[8] * m2.data[8]};
}

static inline mat3 mat3_mul_scalar(mat3 m1, float scalar) {
  return (mat3){.data[0] = m1.data[0] * scalar,
                .data[1] = m1.data[1] * scalar,
                .data[2] = m1.data[2] * scalar,
                .data[3] = m1.data[3] * scalar,
                .data[4] = m1.data[4] * scalar,
                .data[5] = m1.data[5] * scalar,
                .data[6] = m1.data[6] * scalar,
                .data[7] = m1.data[7] * scalar,
                .data[8] = m1.data[8] * scalar};
}

static inline mat3 mat3_add_mat3(mat3 m1, mat3 m2) {
  return (mat3){.data[0] = m1.data[0] + m2.data[0],
                .data[1] = m1.data[1] + m2.data[1],
                .data[2] = m1.data[2] + m2.data[2],
                .data[3] = m1.data[3] + m2.data[3],
                .data[4] = m1.data[4] + m2.data[4],
                .data[5] = m1.data[5] + m2.data[5],
                .data[6] = m1.data[6] + m2.data[6],
                .data[7] = m1.data[7] + m2.data[7],
                .data[8] = m1.data[8] + m2.data[8]};
}

static inline mat3 mat3_orientation(mat3 q1) {
  return (mat3){.data[0] = 1 - (2 * q1.data[2] * q1.data[2] + 2 * q1.data[3] * q1.data[3]),
                .data[1] = 2 * q1.data[1] * q1.data[2] + 2 * q1.data[3] * q1.data[0],
                .data[2] = 2 * q1.data[1] * q1.data[3] - 2 * q1.data[2] * q1.data[0],
                .data[3] = 2 * q1.data[1] * q1.data[2] - 2 * q1.data[3] * q1.data[0],
                .data[4] = 1 - (2 * q1.data[1] * q1.data[1] + 2 * q1.data[3] * q1.data[3]),
                .data[5] = 2 * q1.data[2] * q1.data[3] + 2 * q1.data[1] * q1.data[0],
                .data[6] = 2 * q1.data[1] * q1.data[3] + 2 * q1.data[2] * q1.data[0],
                .data[7] = 2 * q1.data[2] * q1.data[3] - 2 * q1.data[1] * q1.data[0],
                .data[8] = 1 - (2 * q1.data[1] * q1.data[1] + 2 * q1.data[2] * q1.data[2])};
}

static inline mat3 mat3_liner_interpolate(mat3 m1, mat3 m2, float prop) {
  return (mat3){.data[0] = m1.data[0] * (1 - prop) + m2.data[0] * prop,
                .data[1] = m1.data[1] * (1 - prop) + m2.data[1] * prop,
                .data[2] = m1.data[2] * (1 - prop) + m2.data[2] * prop,
                .data[3] = m1.data[3] * (1 - prop) + m2.data[3] * prop,
                .data[4] = m1.data[4] * (1 - prop) + m2.data[4] * prop,
                .data[5] = m1.data[5] * (1 - prop) + m2.data[5] * prop,
                .data[6] = m1.data[6] * (1 - prop) + m2.data[6] * prop,
                .data[7] = m1.data[7] * (1 - prop) + m2.data[7] * prop,
                .data[8] = m1.data[8] * (1 - prop) + m2.data[8] * prop};
}

#endif  // MAT3_H
