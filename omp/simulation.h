
#ifndef __SIMULATION_H
#define __SIMULATION_H

#include "particles.h"

#define CHUNK_SIZE 10

typedef struct {
  particle_t particle_zero;
  long long number_of_collisions;
} simulation_result;

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep, std::vector<particle_t> &par);

#endif // __SIMULATION_H
