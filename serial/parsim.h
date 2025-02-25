
#ifndef __PARSIM_H
#define __PARSIM_H

#include "particles.h"

typedef struct {
  particle_t particle_zero;
  uint64_t number_of_collisions;
} simulation_result;

simulation_result simulation(double side, uint32_t ncside, uint64_t npart,
                             uint64_t nstep,
                             const std::vector<particle_t> &par);

#endif // __PARSIM_H
