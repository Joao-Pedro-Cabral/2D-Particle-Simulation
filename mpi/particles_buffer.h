
#ifndef __PARTICLES_BUFFER_H
#define __PARTICLES_BUFFER_H

#include "particles.h"

typedef struct {
  long long* counters;
  long num_counters;
  long long size;
  long long capacity;
  long long total_num_particles;
  particle_t *particles;
} particles_buffer_t;

void particles_buffer_init(particles_buffer_t *buffer, long num_counters);
void particles_buffer_increment_counter(particles_buffer_t * buffer, long cell);
int particles_buffer_size();
void particles_buffer_realloc(particles_buffer_t* buffer, long long capacity);
void particles_buffer_clean_counters(particles_buffer_t *buffer);
void particles_buffer_clean(particles_buffer_t *buffer);

#endif // __PARTICLES_BUFFER_H
