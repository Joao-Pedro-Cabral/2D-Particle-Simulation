
#ifndef __PARSIM_H
#define __PARSIM_H

#include "particles.h"

typedef struct {
  particle_t particle_zero;
  long long number_of_collisions;
} simulation_result;

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep,
                             std::vector<particle_t> &par);

#endif // __PARSIM_H
