#include <stdio.h>

#include "../include/ubermath/ubermath.h"

int main(int argc, char* argv[]) {
  vec2 pos1 = {.x = 1.0f, .y = 1.0f};
  vec2 pos2 = {.x = 1.5f, .y = 1.5f};
  vec2 pos3 = vec2_add(pos1, pos2);
  vec2 pos4 = pos3;
  pos4 = vec2_add(pos4, pos1);

  return 0;
}