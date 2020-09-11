#include <stdio.h>

#include "../include/ubermath/ubermath.h"

int main(int argc, char* argv[]) {
  vec2 pos1 = {.x = 1.0f, .y = 1.0f};
  vec2 pos2 = {.x = 1.5f, .y = 1.5f};
  vec2 pos3 = vec2_add(pos1, pos2);
  vec2 pos4 = pos3;
  pos4 = vec2_add(pos4, pos1);

  printf("Particle 1 x: %f y: %f\n", pos1.x, pos1.y);
  printf("Particle 2 x: %f y: %f\n", pos2.x, pos2.y);
  printf("Particle 3 x: %f y: %f\n", pos3.x, pos3.y);
  printf("Particle 4 x: %f y: %f\n", pos4.x, pos4.y);

  return 0;
}