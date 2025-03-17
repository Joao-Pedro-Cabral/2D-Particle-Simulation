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

void compute_centers_of_mass(double side, long ncside, long ncside2,
                             cell_t *cells, particle_t *centers) {

  for (long i = 0; i < ncside2; i++) {
    double total_mass = 0.0;
    double weighted_x = 0.0;
    double weighted_y = 0.0;
    long long cell_size = cells[i].size;

    long long iter = (cell_size) - (cell_size & 3);
    long long j;
    __m256d mass_vec = _mm256_setzero_pd();
    __m256d x_vec = _mm256_setzero_pd();
    __m256d y_vec = _mm256_setzero_pd();
    for (j = 0; j < iter; j+=4) {
      __m256d m =  _mm256_loadu_pd(&cells[i].m[j]);
      __m256d x =  _mm256_loadu_pd(&cells[i].x[j]);
      __m256d y =  _mm256_loadu_pd(&cells[i].y[j]);
      mass_vec = _mm256_add_pd(m, mass_vec);
      x_vec = _mm256_fmadd_pd(m, x, x_vec);
      y_vec = _mm256_fmadd_pd(m, y, y_vec);
    }
    total_mass = sum_lanes(mass_vec);
    weighted_x = sum_lanes(x_vec);
    weighted_y = sum_lanes(y_vec);
    while(j < cell_size) {
      total_mass += cells[i].m[j];
      weighted_x += cells[i].m[j] * cells[i].x[j];
      weighted_y += cells[i].m[j] * cells[i].y[j];
      j++;
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
  }
  for (long i = 0; i < ncside + 2; i += ncside + 1) {
    for (long j = 0; j < ncside + 2; j++) {
      long ind = i * (ncside + 2) + j;
      long i2 = (i == 0) ? ncside : 1;
      if (j == 0 || j == ncside + 1) {
        long j2 = (j == 0) ? ncside : 1;
        long ind2 = i2 * (ncside + 2) + j2;
        centers[ind].x =
            (j == 0) ? centers[ind2].x - side : centers[ind2].x + side;
        centers[ind].y =
            (i == 0) ? centers[ind2].y - side : centers[ind2].y + side;
        centers[ind].m = centers[ind2].m;
      } else {
        long ind2 = i2 * (ncside + 2) + j;
        centers[ind].x = centers[ind2].x;
        centers[ind].y =
            (i == 0) ? centers[ind2].y - side : centers[ind2].y + side;
        centers[ind].m = centers[ind2].m;
      }
    }
  }
  for (long j = 0; j < ncside + 2; j += ncside + 1) {
    for (long i = 1; i < ncside + 1; i++) {
      long ind = i * (ncside + 2) + j;
      long j2 = (j == 0) ? ncside : 1;
      long ind2 = i * (ncside + 2) + j2;
      centers[ind].x =
          (j == 0) ? centers[ind2].x - side : centers[ind2].x + side;
      centers[ind].y = centers[ind2].y;
      centers[ind].m = centers[ind2].m;
    }
  }
  debug_centers(ncside, centers);
}

void compute_kinetics(double side, long ncside, long ncside2, cell_t *cells,
                      particle_t *centers) {

  for (long i = 0; i < ncside2; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      __m256d resx_vec = _mm256_setzero_pd();
      __m256d resy_vec = _mm256_setzero_pd();
      double mg = G * cells[i].m[j];
      __m256d mg2 = _mm256_set1_pd(mg);
      __m256d xj =  _mm256_set1_pd(cells[i].x[j]);
      __m256d yj =  _mm256_set1_pd(cells[i].y[j]);
      long long iter = (cell_size) - ((cell_size - (j+1)) & 3);
      long long k;
      for (k = j + 1; k < iter; k+= 4) {
        __m256d xk =  _mm256_loadu_pd(&cells[i].x[k]);
        __m256d yk =  _mm256_loadu_pd(&cells[i].y[k]);
        __m256d dx = _mm256_sub_pd(xk, xj);
        __m256d dy = _mm256_sub_pd(yk, yj);
        __m256d denominator = _mm256_mul_pd(dx, dx);
        denominator = _mm256_fmadd_pd(dy, dy, denominator);
        __m256d sqrt_den = _mm256_sqrt_pd(denominator);
        denominator = _mm256_mul_pd(sqrt_den, denominator);
        __m256d mk = _mm256_loadu_pd(&cells[i].m[k]);
        __m256d numerator = _mm256_mul_pd(mg2, mk);
        __m256d F = _mm256_div_pd(numerator, denominator);
        __m256d forcex = _mm256_mul_pd(dx, F);
        __m256d forcey = _mm256_mul_pd(dy, F);
        __m256d axk =  _mm256_loadu_pd(&cells[i].ax[k]);
        __m256d ayk =  _mm256_loadu_pd(&cells[i].ay[k]);
        axk = _mm256_sub_pd(axk, forcex);
        ayk = _mm256_sub_pd(ayk, forcey);
        _mm256_storeu_pd(&cells[i].ax[k], axk);
        _mm256_storeu_pd(&cells[i].ay[k], ayk);
        resx_vec = _mm256_add_pd(resx_vec, forcex);
        resy_vec = _mm256_add_pd(resy_vec, forcey);
      }
      double resx = sum_lanes(resx_vec);
      double resy = sum_lanes(resy_vec);
      while(k < cell_size) {
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
        k++;
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
                               cell_t *cells) {
  for (long i = 0; i < ncside2; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      long ind = find_cell_c(&cells[i], j, size, ncside);
      if (ind == i)
        continue;
      cell_push_back_c(&cells[ind], &cells[i], j);
      cells[i].ind[j] = -1;
    }
  }
}

void detect_collisions(long ncside2, cell_t *cells) {

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
      __m256i collidedj = _mm256_set1_epi64x(cells[i].collided[j]);
      __m256d epsilon2 = _mm256_set1_pd(EPSILON2);
      __m256d xj =  _mm256_set1_pd(cells[i].x[j]);
      __m256d yj =  _mm256_set1_pd(cells[i].y[j]);
      long long k;
      long long iter = (j) - (j & 3);
      for (k = 0; k < iter; k+= 4) {
        __m256d xk =  _mm256_loadu_pd(&cells[i].x[k]);
        __m256d yk =  _mm256_loadu_pd(&cells[i].y[k]);
        __m256d dx = _mm256_sub_pd(xk, xj);
        __m256d dy = _mm256_sub_pd(yk, yj);
        __m256d distance = _mm256_mul_pd(dx, dx);
        distance = _mm256_fmadd_pd(dy, dy, distance);
        __m256d cmp = _mm256_cmp_pd(distance, epsilon2, _CMP_LT_OQ);
        if(_mm256_testz_pd(cmp, cmp)) {
          continue;
        }
        __m256i near = _mm256_castpd_si256(cmp);
        __m256i collided = _mm256_loadu_si256((__m256i*)&cells[i].collided[k]);
        __m256i mask = _mm256_andnot_si256(collided, near);
        collision += (_mm256_movemask_epi8(mask) != 0);
        _mm256_storeu_si256((__m256i*)&cells[i].collided[k], _mm256_or_si256(near, collided));
      }
      while(k < j) {
        double dx = cells[i].x[j] - cells[i].x[k];
        double dy = cells[i].y[j] - cells[i].y[k];
        double distance = dx * dx + dy * dy;
        if (distance > EPSILON2) {
          k++;
          continue;
        }
        DEBUG("Distance: %.6lf, i: %ld, j: %lld, k: %lld\n", distance, i,
              cells[i].ind[j], cells[i].ind[k]);
        if (cells[i].collided[k] == 0) {
          collision++; 
        }
        cells[i].collided[k] = -1;
        k++;
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
  for (long long i = 0; i < nstep; i++) {
    // printf("--------STEP: %lld --------------\n", i);
    compute_centers_of_mass(side, ncside, ncside2, cells, centers);
    compute_kinetics(side, ncside, ncside2, cells, centers);
    compute_new_particle_cell(size, ncside, ncside2, cells);
    detect_collisions(ncside2, cells);
  }
  res.number_of_collisions = ncollisions;
  res.particle_zero = find_particle_zero(ncside2, cells);
  clean_cells(ncside2, cells);
  free(cells);
  free(centers);
  return res;
}
