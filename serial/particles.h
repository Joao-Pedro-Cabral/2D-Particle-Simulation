
#ifndef __PARTICLES_H
#define __PARTICLES_H

#include <omp.h>
#include <vector>

typedef struct {
  double x;
  double y;
  double vx;
  double vy;
  double m;
} particle_t;

typedef struct cell_t {
  particle_t center;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> m;
  std::vector<long long> ind;
  std::vector<bool> collided;
  omp_lock_t lock;

  void push_back(const particle_t &par, long long i);
  void push_back(const cell_t &cell, long long i);
  void resize(long long size);
  void reserve(long long size);
} cell_t;

typedef struct {
  double x;
  double y;
} vec_t;

#define INDEX(x, y, ncside) ((y + 1) * (ncside + 2) + (x + 1))

long find_particle_cell(const particle_t &par, double size, long ncside);
long find_particle_cell(const cell_t &cell, long long i, double size,
                        long ncside);
void copy_particle(particle_t &par1, const cell_t &cell, long long i);
void copy_particle(cell_t &cell, long long i, long long j);
double calc_squared_distance(const cell_t &cell, long long i, long long j);
void update_position_and_velocity(double side, cell_t &cell, long long i);
vec_t calc_gravitational_force(const cell_t &cell, long long i, long long j);
vec_t calc_gravitational_force(const cell_t &cell, long long i,
                               const particle_t &center);

#endif // __PARTICLES_H
