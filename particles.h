
#ifndef __PARTICLES_H
#define __PARTICLES_H

#include <vector>

typedef struct {
  double x;
  double y;
  double vx;
  double vy;
  double m;
} particle_t;

typedef struct {
  double x;
  double y;
  std::vector<particle_t> par;
} cell_t;

#endif // __PARTICLES_H
