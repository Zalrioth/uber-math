#pragma once
#ifndef UBER_MATH_H
#define UBER_MATH_H

#define UM_PI 3.14159f

#include "ivec3.h"
#include "mat3.h"
#include "mat4.h"
#include "quat.h"
#include "vec2.h"
#include "vec3.h"
#include "vec4.h"

static inline mat4 quaternion_to_mat4(quat rotation);
static inline quat mat4_to_quaternion(mat4 matrix);

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
  dest.vecs[0].data[0] = 1 - 2 * (ySquared + zSquared);
  dest.vecs[0].data[1] = 2 * (xy - zw);
  dest.vecs[0].data[2] = 2 * (xz + yw);
  dest.vecs[0].data[3] = 0;
  dest.vecs[1].data[0] = 2 * (xy + zw);
  dest.vecs[1].data[1] = 1 - 2 * (xSquared + zSquared);
  dest.vecs[1].data[2] = 2 * (yz - xw);
  dest.vecs[1].data[3] = 0;
  dest.vecs[2].data[0] = 2 * (xz - yw);
  dest.vecs[2].data[1] = 2 * (yz + xw);
  dest.vecs[2].data[2] = 1 - 2 * (xSquared + ySquared);
  dest.vecs[2].data[3] = 0;
  dest.vecs[3].data[0] = 0;
  dest.vecs[3].data[1] = 0;
  dest.vecs[3].data[2] = 0;
  dest.vecs[3].data[3] = 1;

  return dest;
}

static inline quat mat4_to_quaternion(mat4 matrix) {
  quat dest;
  float diagonal = matrix.vecs[0].data[0] + matrix.vecs[1].data[1] + matrix.vecs[2].data[2];
  if (diagonal > 0) {
    float w4 = (float)(sqrtf(diagonal + 1.0f) * 2.0f);
    dest.data[3] = w4 / 4.0f;
    dest.data[0] = (matrix.vecs[2].data[1] - matrix.vecs[1].data[2]) / w4;
    dest.data[1] = (matrix.vecs[0].data[2] - matrix.vecs[2].data[0]) / w4;
    dest.data[2] = (matrix.vecs[1].data[0] - matrix.vecs[0].data[1]) / w4;
  } else if ((matrix.vecs[0].data[0] > matrix.vecs[1].data[1]) && (matrix.vecs[0].data[0] > matrix.vecs[2].data[2])) {
    float x4 = (float)(sqrtf(1.0f + matrix.vecs[0].data[0] - matrix.vecs[1].data[1] - matrix.vecs[2].data[2]) * 2.0f);
    dest.data[3] = (matrix.vecs[2].data[1] - matrix.vecs[1].data[2]) / x4;
    dest.data[0] = x4 / 4.0f;
    dest.data[1] = (matrix.vecs[0].data[1] + matrix.vecs[1].data[0]) / x4;
    dest.data[2] = (matrix.vecs[0].data[2] + matrix.vecs[2].data[0]) / x4;
  } else if (matrix.vecs[1].data[1] > matrix.vecs[2].data[2]) {
    float y4 = (float)(sqrtf(1.0f + matrix.vecs[1].data[1] - matrix.vecs[0].data[0] - matrix.vecs[2].data[2]) * 2.0f);
    dest.data[3] = (matrix.vecs[0].data[2] - matrix.vecs[2].data[0]) / y4;
    dest.data[0] = (matrix.vecs[0].data[1] + matrix.vecs[1].data[0]) / y4;
    dest.data[1] = y4 / 4.0f;
    dest.data[2] = (matrix.vecs[1].data[2] + matrix.vecs[2].data[1]) / y4;
  } else {
    float z4 = (float)(sqrtf(1.0f + matrix.vecs[2].data[2] - matrix.vecs[0].data[0] - matrix.vecs[1].data[1]) * 2.0f);
    dest.data[3] = (matrix.vecs[1].data[0] - matrix.vecs[0].data[1]) / z4;
    dest.data[0] = (matrix.vecs[0].data[2] + matrix.vecs[2].data[0]) / z4;
    dest.data[1] = (matrix.vecs[1].data[2] + matrix.vecs[2].data[1]) / z4;
    dest.data[2] = z4 / 4.0f;
  }

  float mag = sqrtf(dest.data[3] * dest.data[3] + dest.data[0] * dest.data[0] + dest.data[1] * dest.data[1] + dest.data[2] * dest.data[2]);
  dest.data[3] /= mag;
  dest.data[0] /= mag;
  dest.data[1] /= mag;
  dest.data[2] /= mag;
}

#endif  // UBER_MATH_H
