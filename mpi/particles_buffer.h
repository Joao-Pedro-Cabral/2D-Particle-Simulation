
#ifndef __PARTICLES_BUFFER_H
#define __PARTICLES_BUFFER_H

#include "cells.h"

typedef struct {
  long long size;
  long long capacity;
  particle_t *particles;
} particles_buffer_t;

void particles_buffer_init(particles_buffer_t *buffer, long long capacity);
void particles_buffer_resize(particles_buffer_t* buffer, long long size);
void particles_buffer_add(particles_buffer_t* buffer, cell_t * cell, long long i);
void particles_buffer_clean(particles_buffer_t *buffer);

#endif // __PARTICLES_BUFFER_H
