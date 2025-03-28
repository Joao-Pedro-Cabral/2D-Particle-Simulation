#include "simulation.h"
#include "cells.h"
#include "debug.h"
#include "init_particles.h"
#include "particles.h"
#include <math.h>
#include <immintrin.h>

static long long ncollisions = 0;

void init_structures(double size, long ncside, long ncside2, long long npart,
                     particle_t *par, cell_t *cells) {
  long long *count = malloc(sizeof(long long) * ncside2);
  for (long i = 0; i < ncside2; i++) {
    count[i] = 0;
  }
  for (long long i = 0; i < npart; i++) {
    count[find_cell_p(&par[i], size, ncside)]++;
  }
  long long min_size = 2*npart / ncside2;
  min_size = (min_size > 10) ? min_size : 10;
  for (long i = 0; i < ncside2; i++) {
    count[i] *= 2;
    long long capacity =
        count[i] > npart ? npart : (count[i] > min_size ? count[i] : min_size);
    cell_init(&cells[i], capacity,
              (i % ncside + 1) + (i / ncside + 1) * (ncside + 2));
  }
  for (long long i = 0; i < npart; i++) {
    cell_push_back_p(&cells[find_cell_p(&par[i], size, ncside)], &par[i], i);
  }
  free(count);
}

void clean_cells(long ncside2, cell_t *cells) {
  for (long i = 0; i < ncside2; i++) {
    cell_clean(&cells[i]);
  }
}

void debug_centers(long ncside, particle_t *centers) {
  for (long i = 0; i < (ncside + 2) * (ncside + 2); i++) {
    DEBUG("Center %ld x: %.6lf y: %.6lf m: %.6lf\n", i, centers[i].x,
          centers[i].y, centers[i].m);
  }
}

static inline double sum_lanes(__m256d vec) {
  vec = _mm256_hadd_pd(vec, vec);
  return ((double*)&vec)[0] + ((double*)&vec)[2];
}

void compute_centers_of_mass(double side, long ncside,
                             cell_t *cells, particle_t *centers, long chunk_size) {
#pragma omp for collapse(2) schedule(dynamic, chunk_size)
  for (long ix = 0; ix < ncside; ix++) {
    for (long iy = 0; iy < ncside; iy++) {
    long i = ix + ncside*iy;
    double total_mass = 0.0;
    double weighted_x = 0.0;
    double weighted_y = 0.0;
    long long cell_size = cells[i].size;

    long long j;
    for (j = 0; j < cell_size; j++) {
      total_mass += cells[i].m[j];
      weighted_x += cells[i].m[j] * cells[i].x[j];
      weighted_y += cells[i].m[j] * cells[i].y[j];
    }

    long ind = cells[i].center;
    centers[ind].m = total_mass;

    if (total_mass > 0) {
      centers[ind].x = weighted_x / total_mass;
      centers[ind].y = weighted_y / total_mass;
    } else {
      centers[ind].x = -2 * side;
      centers[ind].y = -2 * side;
    }

    if(ix == 0 || ix == ncside - 1) {
      long ind2 = (ix == 0) ? ind + ncside : ind - ncside;
      centers[ind2].x = (ix == 0) ? centers[ind].x + side : centers[ind].x - side;
      centers[ind2].y = centers[ind].y;
      centers[ind2].m = centers[ind].m;
      if(iy == 0 || iy == ncside - 1) {
        if(ix == 0) {
          long ind3 = ind + ncside;
          ind3 = (iy == 0) ? ind3 + ncside*(ncside + 2) : ind3 - ncside*(ncside + 2);
          centers[ind3].x = centers[ind].x + side;
          centers[ind3].y = (iy == 0) ? centers[ind].y + side : centers[ind].y - side;
          centers[ind3].m = centers[ind].m;
        }
        if(ix == ncside-1) {
          long ind3 = ind - ncside;
          ind3 = (iy == 0) ? ind3 + ncside*(ncside + 2) : ind3 - ncside*(ncside + 2);
          centers[ind3].x = centers[ind].x - side;
          centers[ind3].y = (iy == 0) ? centers[ind].y + side : centers[ind].y - side;
          centers[ind3].m = centers[ind].m;
        }
      }
    }

    if(iy == 0 || iy == ncside - 1) {
      long ind2 = (iy == 0) ? ind + ncside*(ncside + 2) : ind - ncside*(ncside + 2);
      centers[ind2].x = centers[ind].x;
      centers[ind2].y = (iy == 0) ? centers[ind].y + side : centers[ind].y - side;
      centers[ind2].m = centers[ind].m;
    }
  }
  }
}

void compute_kinetics(double side, long ncside, long ncside2, cell_t *cells,
                      particle_t *centers, long chunk_size) {
#pragma omp for schedule(dynamic, chunk_size)
  for (long i = 0; i < ncside2; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      double mg = G * cells[i].m[j];
      double resx = 0.0;
      double resy = 0.0;
      long long k;
      for (k = j + 1; k < cell_size; k++) {
        double dx = cells[i].x[k] - cells[i].x[j];
        double dy = cells[i].y[k] - cells[i].y[j];
        double denominator = dx * dx + dy * dy;
        denominator *= sqrt(denominator);
        double numerator = mg * cells[i].m[k];
        double F = numerator / denominator;
        double forcex = dx * F;
        double forcey = dy * F;
        cells[i].ax[k] -= forcex;
        cells[i].ay[k] -= forcey;
        resx += forcex;
        resy += forcey;
      }
      cells[i].ax[j] += resx;
      cells[i].ay[j] += resy;
      for (long long k = 0; k < 9; k++) {
        if (k == 4)
          continue;
        long ind = cells[i].center + ((k % 3) - 1) + (k / 3 - 1) * (ncside + 2);
        gravitational_force_pc(&cells[i], j, &centers[ind]);
      }
      cells[i].ax[j] /= cells[i].m[j];
      cells[i].ay[j] /= cells[i].m[j];
      update_position_and_velocity(&cells[i], j, side);
    }
  }
}

void debug_particles(long ncside2, cell_t *cells) {
  for (long i = 0; i < ncside2; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      DEBUG("Particle %lld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: %.6f, "
            "ax: %.6f, ay: %.6f, collided: %lld, cell: %ld\n",
            cells[i].ind[j], cells[i].m[j], cells[i].x[j], cells[i].y[j],
            cells[i].vx[j], cells[i].vy[j], cells[i].ax[j], cells[i].ay[j],
            cells[i].collided[j], i);
    }
  }
}

void compute_new_particle_cell(double size, long ncside, long ncside2,
                               cell_t *cells, long chunk_size) {
#pragma omp for schedule(dynamic, chunk_size)
  for (long i = 0; i < ncside2; i++) {
    omp_set_lock(&cells[i].lock);
    long long cell_size = cells[i].size;
    omp_unset_lock(&cells[i].lock);
    for (long long j = 0; j < cell_size; j++) {
      long ind = find_cell_c(&cells[i], j, size, ncside);
      if (ind == i)
        continue;
      omp_set_lock(&cells[ind].lock);
      cell_push_back_c(&cells[ind], &cells[i], j);
      omp_unset_lock(&cells[ind].lock);
      cells[i].ind[j] = -1;
    }
  }
}

void detect_collisions(long ncside2, cell_t *cells, long chunk_size) {
#pragma omp for reduction(+ : ncollisions) schedule(dynamic, chunk_size)
  for (long i = 0; i < ncside2; i++) {
    long long cell_collisions = 0;
    long long cell_size = cells[i].size;
    for (long long j = cell_size - 1; j >= 0; j--) {
      if (cells[i].ind[j] == -1) {
        DEBUG("Removing mass null Particle %lld from cells[%ld].par[%lld]\n",
              cells[i].ind[j], i, j);
        cell_size--;
        cell_remove_particle(&cells[i], j);
      }
    }
    for (long long j = cell_size - 1; j >= 0; j--) {
      long long collision = 0;
      long long k;
      for (k = 0; k < j; k++) {
        double dx = cells[i].x[j] - cells[i].x[k];
        double dy = cells[i].y[j] - cells[i].y[k];
        double distance = dx * dx + dy * dy;
        if (distance > EPSILON2) {
          continue;
        }
        DEBUG("Distance: %.6lf, i: %ld, j: %lld, k: %lld\n", distance, i,
              cells[i].ind[j], cells[i].ind[k]);
        if (cells[i].collided[k] == 0) {
          collision++; 
        }
        cells[i].collided[k] = -1;
      }
      if(collision > 0 && (cells[i].collided[j] == 0)) {
        cells[i].collided[j] = -1;
        cell_collisions ++;
      }
      if (cells[i].collided[j] == -1) {
        DEBUG("Removing Particle %lld from cells[%ld].par[%lld]\n",
              cells[i].ind[j], i, j);
        cell_size--;
        cell_remove_particle(&cells[i], j);
      }
    }
    ncollisions += cell_collisions;
    cell_resize(&cells[i], cell_size);
  }
  DEBUG("Number of collisions %lld\n", ncollisions);
  debug_particles(ncside2, cells);
}

particle_t find_particle_zero(long ncside2, cell_t *cells) {
  particle_t par0;
  for (long i = 0; i < ncside2; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      if (cells[i].ind[j] == 0) {
        par0.x = cells[i].x[j];
        par0.y = cells[i].y[j];
        par0.vx = cells[i].vx[j];
        par0.vy = cells[i].vy[j];
        par0.m = cells[i].m[j];
        return par0;
      }
    }
  }
  ERROR("Particle 0 not found\n");
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep, particle_t *par) {
  double size = side / ncside;
  long ncside2 = ncside * ncside;
  cell_t *cells = malloc(sizeof(cell_t) * ncside2);
  particle_t *centers =
      malloc(sizeof(particle_t) * (ncside + 2) * (ncside + 2));
  simulation_result res;
  init_structures(size, ncside, ncside2, npart, par, cells);
  long chunk_size = (ncside2 >= 4*omp_get_max_threads()) ? 2 : 1;
#pragma omp parallel
  for (long long i = 0; i < nstep; i++) {
    // printf("--------STEP: %lld --------------\n", i);
    compute_centers_of_mass(side, ncside, cells, centers, chunk_size);
    compute_kinetics(side, ncside, ncside2, cells, centers, chunk_size);
    compute_new_particle_cell(size, ncside, ncside2, cells, chunk_size);
    detect_collisions(ncside2, cells, chunk_size);
  }
  res.number_of_collisions = ncollisions;
  res.particle_zero = find_particle_zero(ncside2, cells);
  clean_cells(ncside2, cells);
  free(cells);
  free(centers);
  return res;
}
