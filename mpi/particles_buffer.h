
#ifndef __PARTICLES_BUFFER_H
#define __PARTICLES_BUFFER_H

#include "particles.h"

typedef struct {
  particle_t *par;
  long long* size;
  long long capacity;
} particles_buffer_t;

void particles_buffer_init(particles_buffer_t *buffer, long long capacity);
void particles_buffer_resize(particles_buffer_t *buffer, long long size);
void cell_remove_particle(cell_t *cell, long long i);
void particles_buffer_clean(particles_buffer_t *buffer);

#endif // __PARTICLES_BUFFER_H
