#pragma once
#ifndef MAT4_H
#define MAT4_H

#include "ubermathcommon.h"

struct mat4 {
  union {
    struct
    {
      float m11, m21, m31, m41,
          m12, m22, m32, m42,
          m13, m23, m33, m43,
          m14, m24, m34, m44;
    };
    float data[16];
    alignas(16) __m128 sse_data[4];
  };
};

#endif  // MAT4_H
