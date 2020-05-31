#include <stdio.h>

#include "../include/ubermath/ubermath.h"

struct Particles {
  size_t size;
  size_t back;
  vec2_soa vec2_entities;
};

void particles_init(struct Particles* particles, size_t size) {
  particles->size = size;
  particles->back = size - 1;
  vec2_soa_init(&particles->vec2_entities);
  vec2_soa_resize(&particles->vec2_entities, 0, size);

  // Set position of each particle
  for (int particle_num = 0; particle_num < size; particle_num++) {
    particles->vec2_entities.x[particle_num] = (float)rand() / RAND_MAX;
    particles->vec2_entities.y[particle_num] = (float)rand() / RAND_MAX;
    printf("Particle %d position is x: %f y: %f\n", particle_num, particles->vec2_entities.x[particle_num], particles->vec2_entities.y[particle_num]);
  }
  printf("------------------------------------------------------------------------------------------\n");
}

void particles_delete(struct Particles* particles) {
  vec2_soa_delete(&particles->vec2_entities);
}

void particles_remove(struct Particles* particles, size_t position) {
  particles->vec2_entities.x[position] = particles->vec2_entities.x[particles->back];
  particles->vec2_entities.x[position] = particles->vec2_entities.x[particles->back];
  particles->back--;
  vec2_soa_resize(&particles->vec2_entities, particles->size, particles->size - 1);
  particles->size--;
}

// Apply gravity force in groups of max simd support sse4/neon is 4 avx2 is 8
void particles_apply_gravity(struct Particles* particles) {
  simd_float_max gravity = simd_float_max_set1(0.5f);
  for (int iter_num = 0; iter_num < vec2_soa_iterations(particles->size); iter_num++)
    particles->vec2_entities.simd_data_y[iter_num] = vec2_soa_add(particles->vec2_entities.simd_data_y[iter_num], gravity);

  for (int particle_num = 0; particle_num < particles->size; particle_num++)
    printf("Particle %d position after gravity addition is x: %f y: %f\n", particle_num, particles->vec2_entities.x[particle_num], particles->vec2_entities.y[particle_num]);
  printf("------------------------------------------------------------------------------------------\n");
}

int main(int argc, char* argv[]) {
  size_t total_particles = 20;

  struct Particles particles = {0};
  particles_init(&particles, total_particles);

  particles_apply_gravity(&particles);

  // Remove half the particles
  size_t half_particles = particles.size / 2;
  for (int loop_num = 0; loop_num < half_particles; loop_num++)
    particles_remove(&particles, rand() % (particles.size));

  particles_apply_gravity(&particles);

  particles_delete(&particles);

  return 0;
}