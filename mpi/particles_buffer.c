
#include "particles_buffer.h"
#include <stdlib.h>

void particles_buffer_init(particles_buffer_t *buffer, long long capacity) {
  buffer->particles = (particle_t *)malloc(capacity * sizeof(particle_t));
  buffer->size = 0;
  buffer->capacity = capacity;
}

void particles_buffer_resize(particles_buffer_t *buffer, long long size) {
  buffer->size = size;
  if (size > buffer->capacity) {
    buffer->capacity = 2 * size;
    buffer->particles = (particle_t *)realloc(
        buffer->particles, buffer->capacity * sizeof(particle_t));
  }
}

void particles_buffer_add_p(particles_buffer_t *buffer, double x,
        double y, double vx, double vy, double m, long long ind) {
  buffer->particles[buffer->size].x = x;
  buffer->particles[buffer->size].y = y;
  buffer->particles[buffer->size].vx = vx;
  buffer->particles[buffer->size].vy = vy;
  buffer->particles[buffer->size].m = m;
  buffer->particles[buffer->size].ind = ind;
  buffer->size++;
  particles_buffer_resize(buffer, buffer->size);
}

void particles_buffer_add_c(particles_buffer_t *buffer, cell_t *cell,
                          long long i) {
  buffer->particles[buffer->size].x = cell->x[i];
  buffer->particles[buffer->size].y = cell->y[i];
  buffer->particles[buffer->size].vx = cell->vx[i];
  buffer->particles[buffer->size].vy = cell->vy[i];
  buffer->particles[buffer->size].m = cell->m[i];
  buffer->particles[buffer->size].ind = cell->ind[i];
  buffer->size++;
  particles_buffer_resize(buffer, buffer->size);
}

void particles_buffer_clean(particles_buffer_t *buffer) {
  free(buffer->particles);
  buffer->particles = NULL;
  buffer->size = 0;
  buffer->capacity = 0;
}
