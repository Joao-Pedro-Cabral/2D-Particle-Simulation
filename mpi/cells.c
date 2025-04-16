
#include "cells.h"
#include "init_particles.h"
#include <math.h>
#include <stdlib.h>

void cell_init(cell_t *cell, long long capacity, long center) {
  cell->x = (double *)malloc(capacity * sizeof(double));
  cell->y = (double *)malloc(capacity * sizeof(double));
  cell->vx = (double *)malloc(capacity * sizeof(double));
  cell->vy = (double *)malloc(capacity * sizeof(double));
  cell->ax = (double *)malloc(capacity * sizeof(double));
  cell->ay = (double *)malloc(capacity * sizeof(double));
  cell->m = (double *)malloc(capacity * sizeof(double));
  cell->ind = (long long *)malloc(capacity * sizeof(long long));
  cell->collided = (long long *)malloc(capacity * sizeof(long long));
  cell->size = 0;
  cell->capacity = capacity;
  cell->center = center;
  omp_init_lock(&cell->lock);
}

void cell_resize(cell_t *cell, long long size) {
  cell->size = size;
  if (size >= cell->capacity) {
    long long capacity = (cell->capacity == 0) ? 1 : cell->capacity * 2;
    cell->x = (double *)realloc(cell->x, capacity * sizeof(double));
    cell->y = (double *)realloc(cell->y, capacity * sizeof(double));
    cell->vx = (double *)realloc(cell->vx, capacity * sizeof(double));
    cell->vy = (double *)realloc(cell->vy, capacity * sizeof(double));
    cell->ax = (double *)realloc(cell->ax, capacity * sizeof(double));
    cell->ay = (double *)realloc(cell->ay, capacity * sizeof(double));
    cell->m = (double *)realloc(cell->m, capacity * sizeof(double));
    cell->ind = (long long *)realloc(cell->ind, capacity * sizeof(long long));
    cell->collided =
        (long long *)realloc(cell->collided, capacity * sizeof(long long));
    cell->capacity = capacity;
  }
}

void cell_push_back_c(cell_t *cell, cell_t *cell2, long long i) {
  cell->x[cell->size] = cell2->x[i];
  cell->y[cell->size] = cell2->y[i];
  cell->vx[cell->size] = cell2->vx[i];
  cell->vy[cell->size] = cell2->vy[i];
  cell->ax[cell->size] = cell2->ax[i];
  cell->ay[cell->size] = cell2->ay[i];
  cell->m[cell->size] = cell2->m[i];
  cell->ind[cell->size] = cell2->ind[i];
  cell->collided[cell->size] = cell2->collided[i];
  cell->size++;
  cell_resize(cell, cell->size);
}

void cell_push_back_p(cell_t *cell, particle_t *par) {
  cell->x[cell->size] = par->x;
  cell->y[cell->size] = par->y;
  cell->vx[cell->size] = par->vx;
  cell->vy[cell->size] = par->vy;
  cell->ax[cell->size] = 0.0;
  cell->ay[cell->size] = 0.0;
  cell->m[cell->size] = par->m;
  cell->ind[cell->size] = par->ind;
  cell->collided[cell->size] = 0;
  cell->size++;
  cell_resize(cell, cell->size);
}

void cell_remove_particle(cell_t *cell, long long i) {
  cell->size = cell->size - 1;
  cell->x[i] = cell->x[cell->size];
  cell->y[i] = cell->y[cell->size];
  cell->vx[i] = cell->vx[cell->size];
  cell->vy[i] = cell->vy[cell->size];
  cell->ax[i] = cell->ax[cell->size];
  cell->ay[i] = cell->ay[cell->size];
  cell->m[i] = cell->m[cell->size];
  cell->ind[i] = cell->ind[cell->size];
  cell->collided[i] = cell->collided[cell->size];
}

void cell_clean(cell_t *cell) {
  free(cell->x);
  free(cell->y);
  free(cell->vx);
  free(cell->vy);
  free(cell->ax);
  free(cell->ay);
  free(cell->m);
  free(cell->ind);
  free(cell->collided);
  cell->x = NULL;
  cell->y = NULL;
  cell->vx = NULL;
  cell->vy = NULL;
  cell->ax = NULL;
  cell->ay = NULL;
  cell->m = NULL;
  cell->ind = NULL;
  cell->collided = NULL;
  cell->size = 0;
  cell->capacity = 0;
  cell->center = 0;
  omp_destroy_lock(&cell->lock);
}

double squared_distance(cell_t *cell, long long i, long long j) {
  double dx = cell->x[i] - cell->x[j];
  double dy = cell->y[i] - cell->y[j];
  return dx * dx + dy * dy;
}

void update_position_and_velocity(cell_t *cell, long long i, double side) {
  cell->x[i] += DELTAT * cell->vx[i] + 0.5 * (DELTAT * DELTAT) * cell->ax[i];
  if (cell->x[i] > side)
    cell->x[i] -= side;
  if (cell->x[i] < 0)
    cell->x[i] += side;
  cell->y[i] += DELTAT * cell->vy[i] + 0.5 * (DELTAT * DELTAT) * cell->ay[i];
  if (cell->y[i] > side)
    cell->y[i] -= side;
  if (cell->y[i] < 0)
    cell->y[i] += side;
  cell->vx[i] += DELTAT * cell->ax[i];
  cell->vy[i] += DELTAT * cell->ay[i];
  cell->ax[i] = 0.0;
  cell->ay[i] = 0.0;
}

void gravitational_force_pp(cell_t *cell, long long i, long long j) {
  double dx = cell->x[j] - cell->x[i];
  double dy = cell->y[j] - cell->y[i];
  double denominator = dx * dx + dy * dy;
  denominator *= sqrt(denominator);
  double numerator = G * cell->m[i] * cell->m[j];
  double F = numerator / denominator;
  double force_x = dx * F;
  double force_y = dy * F;
  cell->ax[i] += force_x;
  cell->ay[i] += force_y;
  cell->ax[j] -= force_x;
  cell->ay[j] -= force_y;
}

void gravitational_force_pc(cell_t *cell, long long i, center_t *center) {
  double dx = center->x - cell->x[i];
  double dy = center->y - cell->y[i];
  double denominator = dx * dx + dy * dy;
  denominator *= sqrt(denominator);
  double numerator = G * cell->m[i] * center->m;
  double F = numerator / denominator;
  double force_x = dx * F;
  double force_y = dy * F;
  cell->ax[i] += force_x;
  cell->ay[i] += force_y;
}

void copy_center(center_t *center1, const center_t *center2) {
  center1->x = center2->x;
  center1->y = center2->y;
  center1->m = center2->m;
}
