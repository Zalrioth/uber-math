#pragma once
#ifndef QUAT_H
#define QUAT_H

#include <float.h>

#include "ubermathcommon.h"
#include "vec3.h"
#include "vec4.h"

#define QUAT_ZERO \
  (quat) { .data[0] = 0.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f }

#define QUAT_DEFAULT \
  (quat) { .data[0] = 1.0f, .data[1] = 0.0f, .data[2] = 0.0f, .data[3] = 0.0f }

typedef struct quat {
  union {
    struct {
      float x, y, z, w;
    };
    struct {
      float r, i, j, k;
    };
    //vec4 vec;
    float data[4];
    //simd_align_max float data[4];
  };
} quat;

static inline quat quaternion_set(float r, float i, float j, float k);
static inline quat quaternion_normalise(quat q1);
static inline quat quaternion_mul(quat q1, quat q2);
static inline quat quaternion_add_scaled_vector(quat q1, vec3 v1, float scale);
static inline quat quaternion_rotate_by_vector(quat q1, vec3 v1);
static inline quat quat_interpolate_linear(quat a, quat b, float blend);

static inline quat quaternion_set(float r, float i, float j, float k) {
  return (quat){.data[0] = r, .data[1] = i, .data[2] = j, .data[3] = k};
}

static inline quat quaternion_normalise(quat q1) {
  float d = q1.data[0] * q1.data[0] + q1.data[1] * q1.data[1] + q1.data[2] * q1.data[2] + q1.data[3] * q1.data[3];

  if (d < FLT_EPSILON)
    return (quat){.data[0] = 1, .data[1] = q1.data[1], .data[2] = q1.data[2], .data[3] = q1.data[3]};

  d = ((float)1.0) / sqrtf(d);
  return (quat){.data[0] = q1.data[0] * d, .data[1] = q1.data[1] * d, .data[2] = q1.data[2] * d, .data[3] = q1.data[3] * d};
}

static inline quat quaternion_mul(quat q1, quat q2) {
  return (quat){.data[0] = q1.data[0] * q2.data[0] - q1.data[1] * q2.data[1] - q1.data[2] * q2.data[2] - q1.data[3] * q2.data[3],
                .data[1] = q1.data[0] * q2.data[1] + q1.data[1] * q2.data[0] + q1.data[2] * q2.data[3] - q1.data[3] * q2.data[2],
                .data[2] = q1.data[0] * q2.data[2] + q1.data[2] * q2.data[0] + q1.data[3] * q2.data[1] - q1.data[1] * q2.data[3],
                .data[3] = q1.data[0] * q2.data[3] + q1.data[3] * q2.data[0] + q1.data[1] * q2.data[2] - q1.data[2] * q2.data[1]};
}

static inline quat quaternion_add_scaled_vector(quat q1, vec3 v1, float scale) {
  quat temp = quaternion_mul(q1, (quat){.data[0] = 0, .data[1] = v1.data[0] * scale, .data[2] = v1.data[1] * scale, .data[3] = v1.data[2] * scale});
  return (quat){.data[0] = q1.r + temp.r * 0.5f, .data[1] = q1.i + temp.i * 0.5f, .data[2] = q1.j + temp.j * 0.5f, .data[3] = q1.k + temp.k * 0.5f};
}

static inline quat quaternion_rotate_by_vector(quat q1, vec3 v1) {
  return quaternion_mul(q1, (quat){.data[0] = 0, .data[1] = v1.data[0], .data[2] = v1.data[1], .data[3] = v1.data[2]});
}

static inline quat quat_interpolate_linear(quat a, quat b, float blend) {
  quat dest;
  float dot = a.data[3] * b.data[3] + a.data[0] * b.data[0] + a.data[1] * b.data[1] + a.data[2] * b.data[2];
  float blendI = 1.0f - blend;
  if (dot < 0) {
    dest.data[3] = blendI * a.data[3] + blend * -b.data[3];
    dest.data[0] = blendI * a.data[0] + blend * -b.data[0];
    dest.data[1] = blendI * a.data[1] + blend * -b.data[1];
    dest.data[2] = blendI * a.data[2] + blend * -b.data[2];
  } else {
    dest.data[3] = blendI * a.data[3] + blend * b.data[3];
    dest.data[0] = blendI * a.data[0] + blend * b.data[0];
    dest.data[1] = blendI * a.data[1] + blend * b.data[1];
    dest.data[2] = blendI * a.data[2] + blend * b.data[2];
  }

  float mag = sqrtf(dest.data[3] * dest.data[3] + dest.data[0] * dest.data[0] + dest.data[1] * dest.data[1] + dest.data[2] * dest.data[2]);
  dest.data[3] /= mag;
  dest.data[0] /= mag;
  dest.data[1] /= mag;
  dest.data[2] /= mag;

  return dest;
}

#endif  // QUAT_H
