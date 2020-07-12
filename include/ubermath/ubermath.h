#pragma once
#ifndef UBER_MATH_H
#define UBER_MATH_H

//#define UM_PI 3.14159f
#define UM_PI (float)3.14159265358979323846264338327950288

#include "ivec3.h"
#include "mat3.h"
#include "mat4.h"
#include "quat.h"
#include "vec2.h"
#include "vec3.h"
#include "vec4.h"

static inline float degree_to_radian(float degree);
static inline float radian_to_degree(float radian);
static inline vec4 vec3_to_vec4(vec3 v1);
static inline mat4 quaternion_to_mat4(quat rotation);
static inline quat mat4_to_quaternion(mat4 matrix);

static inline float degree_to_radian(float degree) {
  return degree * UM_PI / 180.0f;
}

static inline float radian_to_degree(float radian) {
  return radian * 180.0f / UM_PI;
}

static inline vec4 vec3_to_vec4(vec3 v1) {
  return (vec4){.data[0] = v1.data[0], .data[1] = v1.data[1], .data[2] = v1.data[2], 0.0f};
}

static inline mat4 quaternion_to_mat4(quat rotation) {
  float xy = rotation.data[0] * rotation.data[1];
  float xz = rotation.data[0] * rotation.data[2];
  float xw = rotation.data[0] * rotation.data[3];
  float yz = rotation.data[1] * rotation.data[2];
  float yw = rotation.data[1] * rotation.data[3];
  float zw = rotation.data[2] * rotation.data[3];
  float xSquared = rotation.data[0] * rotation.data[0];
  float ySquared = rotation.data[1] * rotation.data[1];
  float zSquared = rotation.data[2] * rotation.data[2];

  mat4 dest;
  dest.m00 = 1 - 2 * (ySquared + zSquared);
  dest.m01 = 2 * (xy - zw);
  dest.m02 = 2 * (xz + yw);
  dest.m03 = 0;
  dest.m10 = 2 * (xy + zw);
  dest.m11 = 1 - 2 * (xSquared + zSquared);
  dest.m12 = 2 * (yz - xw);
  dest.m13 = 0;
  dest.m20 = 2 * (xz - yw);
  dest.m21 = 2 * (yz + xw);
  dest.m22 = 1 - 2 * (xSquared + ySquared);
  dest.m23 = 0;
  dest.m30 = 0;
  dest.m31 = 0;
  dest.m32 = 0;
  dest.m33 = 1;

  return dest;
}

static inline quat mat4_to_quaternion(mat4 matrix) {
  quat dest;
  float diagonal = matrix.m00 + matrix.m11 + matrix.m22;
  if (diagonal > 0) {
    float w4 = (float)(sqrtf(diagonal + 1.0f) * 2.0f);
    dest.data[3] = w4 / 4.0f;
    dest.data[0] = (matrix.m21 - matrix.m12) / w4;
    dest.data[1] = (matrix.m02 - matrix.m20) / w4;
    dest.data[2] = (matrix.m10 - matrix.m01) / w4;
  } else if ((matrix.m00 > matrix.m11) && (matrix.m00 > matrix.m22)) {
    float x4 = (float)(sqrtf(1.0f + matrix.m00 - matrix.m11 - matrix.m22) * 2.0f);
    dest.data[3] = (matrix.m21 - matrix.m12) / x4;
    dest.data[0] = x4 / 4.0f;
    dest.data[1] = (matrix.m01 + matrix.m10) / x4;
    dest.data[2] = (matrix.m02 + matrix.m20) / x4;
  } else if (matrix.m11 > matrix.m22) {
    float y4 = (float)(sqrtf(1.0f + matrix.m11 - matrix.m00 - matrix.m22) * 2.0f);
    dest.data[3] = (matrix.m02 - matrix.m20) / y4;
    dest.data[0] = (matrix.m01 + matrix.m10) / y4;
    dest.data[1] = y4 / 4.0f;
    dest.data[2] = (matrix.m12 + matrix.m21) / y4;
  } else {
    float z4 = (float)(sqrtf(1.0f + matrix.m22 - matrix.m00 - matrix.m11) * 2.0f);
    dest.data[3] = (matrix.m10 - matrix.m01) / z4;
    dest.data[0] = (matrix.m02 + matrix.m20) / z4;
    dest.data[1] = (matrix.m12 + matrix.m21) / z4;
    dest.data[2] = z4 / 4.0f;
  }

  float mag = sqrtf(dest.data[3] * dest.data[3] + dest.data[0] * dest.data[0] + dest.data[1] * dest.data[1] + dest.data[2] * dest.data[2]);
  dest.data[3] /= mag;
  dest.data[0] /= mag;
  dest.data[1] /= mag;
  dest.data[2] /= mag;

  return dest;
}

#endif  // UBER_MATH_H
