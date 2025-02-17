
#ifndef __INIT_PARTICLES_H
#define __INIT_PARTICLES_H

#define _USE_MATH_DEFINES
#define G 6.67408e-11
#define EPSILON2 (0.005*0.005)
#define DELTAT 0.1

#include "particles.h"

void init_r4uni(int input_seed);

double rnd_uniform01();

double rnd_normal01();

void init_particles(long seed, double side, long ncside, long long n_part, particle_t *par);

#endif // __INIT_PARTICLES_H
