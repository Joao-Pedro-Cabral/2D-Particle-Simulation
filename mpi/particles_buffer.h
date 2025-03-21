
#ifndef __PARTICLES_BUFFER_H
#define __PARTICLES_BUFFER_H

#include "particles.h"

typedef struct {
  particle_t *particles;
  long long* counters;
  long num_counters;
  long long capacity;
} particles_buffer_t;

void particles_buffer_init(particles_buffer_t *buffer, long num_counters, long long capacity);
void particles_buffer_increment_counter(particles_buffer_t * buffer, long cell);
void particles_buffer_clean_counters(particles_buffer_t *buffer);
void particles_buffer_add(particles_buffer_t* buffer, long cell);
void particles_buffer_clean(particles_buffer_t *buffer);

#endif // __PARTICLES_BUFFER_H
