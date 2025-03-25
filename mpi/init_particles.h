
#ifndef __INIT_PARTICLES_H
#define __INIT_PARTICLES_H

#define _USE_MATH_DEFINES
#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

#include "particles.h"
#include "particles_buffer.h"
#include "communication_buffers.h"

void init_r4uni(int input_seed);

double rnd_uniform01();

double rnd_normal01();

void init_particles(long seed, double side, long ncside, int id, long long n_part,
                    particles_buffer_t *par, communication_buffers_t * buffers);

#endif // __INIT_PARTICLES_H
