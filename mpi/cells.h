
#ifndef __CELL_H
#define __CELL_H

#include "particles.h"
#include <omp.h>

typedef struct {
  double x;
  double y;
  double m;
} center_t;

typedef struct {
  long center;
  double *x;
  double *y;
  double *vx;
  double *vy;
  double *ax;
  double *ay;
  double *m;
  long long *ind;
  long long *collided;
  long long size;
  long long capacity;
} cell_t;

void cell_init(cell_t *cell, long long capacity, long center);
void cell_resize(cell_t *cell, long long capacity);
void cell_push_back_c(cell_t *cell, cell_t *cell2, long long i);
void cell_push_back_p(cell_t *cell, particle_t *par);
void cell_remove_particle(cell_t *cell, long long i);
void cell_clean(cell_t *cell);
double squared_distance(cell_t *cell, long long i, long long j);
void update_position_and_velocity(cell_t *cell, long long i, double side);
void gravitational_force_pp(cell_t *cell, long long i, long long j);
void gravitational_force_pc(cell_t *cell, long long i, center_t *center);
void copy_center(center_t * center1, const center_t * center2);

#endif // __CELL_H
