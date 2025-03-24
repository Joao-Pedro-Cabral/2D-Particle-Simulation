
#ifndef __SIMULATION_H
#define __SIMULATION_H

#include "particles.h"
#include "particles_buffer.h"

typedef struct {
  particle_t particle_zero;
  long long number_of_collisions;
} simulation_result;

simulation_result simulation(double side, long ncside, long long npart, int id,
                             int p, long long nstep, particles_buffer_t *par);

#endif // __SIMULATION_H
