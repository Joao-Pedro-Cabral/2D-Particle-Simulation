
#ifndef __PARTICLES_H
#define __PARTICLES_H

#include <vector>
#include <omp.h>

typedef struct {
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  double m;
  long long ind;
  bool collided;
} particle_t;

typedef struct {
  particle_t center;
  std::vector<particle_t> par;
  omp_lock_t lock;
} cell_t;

typedef struct {
  double x;
  double y;
} vec_t;

#define INDEX(x, y, ncside) ((y + 1) * (ncside + 2) + (x + 1))

long find_particle_cell(const particle_t &par, double size, long ncside);
void copy_particle(particle_t &par1, const particle_t &par2);
void remove_and_swap(std::vector<cell_t> &cells, long i, long long j);
double calc_squared_distance(const particle_t &par1, const particle_t &par2);
void update_position_and_velocity(double side, particle_t &par);
vec_t calc_gravitational_force(const particle_t &par1, const particle_t &par2);

#endif // __PARTICLES_H
