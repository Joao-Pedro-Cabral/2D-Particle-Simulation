
#ifndef __PARTICLES_H
#define __PARTICLES_H

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
  long long n_part;
  long long c_part;
  particle_t *par;
} cell_t;

#endif // __PARTICLES_H
