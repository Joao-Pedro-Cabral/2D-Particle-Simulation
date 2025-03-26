#include "simulation.h"
#include "cells.h"
#include "communication_buffers.h"
#include "debug.h"
#include "init_particles.h"
#include "nodes.h"
#include "particles.h"
#include <immintrin.h>
#include <math.h>
#include <mpi.h>

static long long ncollisions = 0;

void init_block(double size, long ncside, long long npart, int id,
                 particles_buffer_t *par, cell_t *cells,
                 communication_buffers_t *buffers) {
  long long *count = malloc(sizeof(long long) * buffers->size);
  for (long i = 0; i < buffers->size; i++) {
    count[i] = 0;
  }
  for (long long i = 0; i < par->size; i++) {
    long cell = find_cell_p(buffers, &par->particles[i], size);
    count[cell]++;
  }
  long long min_size = 2 * npart / (ncside * ncside);
  for (long i = 0; i < buffers->size; i++) {
    count[i] *= 2;
    long long capacity =
        count[i] > npart ? npart : (count[i] > min_size ? count[i] : min_size);
    cell_init(&cells[i], capacity,
              ((i % buffers->lens[1]) + 1) + ((i / buffers->lens[1]) + 1) * (buffers->lens[1] + 2));
  }
  for (long long i = 0; i < par->size; i++) {
    long cell = find_cell_p(buffers, &par->particles[i], size);
    cell_push_back_p(&cells[cell], &par->particles[i]);
  }
  free(count);
  particles_buffer_clean(par);
  communication_buffers_init(buffers, ncside, id, npart);
}

void clean_blocks(cell_t *cells, center_t *centers, communication_buffers_t *buffers) {
  for (long i = 0; i < buffers->size; i++) {
    cell_clean(&cells[i]);
  }
  free(cells);
  free(centers);
  communication_buffers_clean(buffers);
}

void debug_centers(int id, long ncenters, center_t *centers) {
  for (long i = 0; i < ncenters; i++) {
    DEBUG("From %d: Center %ld x: %.6lf y: %.6lf m: %.6lf\n", id, i,
          centers[i].x, centers[i].y, centers[i].m);
  }
}

void debug_particles(int id, long ncells, cell_t *cells) {
  for (long i = 0; i < ncells; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      DEBUG("From %d: Particle %lld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: "
            "%.6f, "
            "ax: %.6f, ay: %.6f, collided: %lld, cell: %ld\n",
            id, cells[i].ind[j], cells[i].m[j], cells[i].x[j], cells[i].y[j],
            cells[i].vx[j], cells[i].vy[j], cells[i].ax[j], cells[i].ay[j],
            cells[i].collided[j], cells[i].center);
    }
  }
}

static inline double sum_lanes(__m256d vec) {
  vec = _mm256_hadd_pd(vec, vec);
  return ((double *)&vec)[0] + ((double *)&vec)[2];
}

void compute_centers_of_mass(double side, int id,
                             cell_t *cells, center_t *centers,
                             communication_buffers_t *buffers) {

  for (long iy = 0; iy < buffers->lens[0]; iy++) {
    for (long ix = 0; ix < buffers->lens[1]; ix++) {
      long i = ix + buffers->lens[1] * iy;
      double total_mass = 0.0;
      double weighted_x = 0.0;
      double weighted_y = 0.0;
      long long cell_size = cells[i].size;

      long long iter = (cell_size) - (cell_size & 3);
      long long j;
      __m256d mass_vec = _mm256_setzero_pd();
      __m256d x_vec = _mm256_setzero_pd();
      __m256d y_vec = _mm256_setzero_pd();
      for (j = 0; j < iter; j += 4) {
        __m256d m = _mm256_loadu_pd(&cells[i].m[j]);
        __m256d x = _mm256_loadu_pd(&cells[i].x[j]);
        __m256d y = _mm256_loadu_pd(&cells[i].y[j]);
        mass_vec = _mm256_add_pd(m, mass_vec);
        x_vec = _mm256_fmadd_pd(m, x, x_vec);
        y_vec = _mm256_fmadd_pd(m, y, y_vec);
      }
      total_mass = sum_lanes(mass_vec);
      weighted_x = sum_lanes(x_vec);
      weighted_y = sum_lanes(y_vec);
      while (j < cell_size) {
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

      if (ix == 0) {
        copy_center(&buffers->send_centers[3][iy], &centers[ind]);
        if(iy == 0) {
          copy_center(&buffers->send_centers[0][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 0]);
        }
        if (iy == buffers->lens[0] - 1) {
          copy_center(&buffers->send_centers[5][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 5]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 3]);
        }
      }
      if (ix == buffers->lens[1] - 1) {
        copy_center(&buffers->send_centers[4][iy], &centers[ind]);
        if(iy == 0) {
          copy_center(&buffers->send_centers[2][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 2]);
        }
        if (iy == buffers->lens[0] - 1) {
          copy_center(&buffers->send_centers[7][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 7]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 4]);
        }
      }
      if(iy == 0) {
        copy_center(&buffers->send_centers[1][ix], &centers[ind]);
        if(ix == buffers->lens[1] - 1) {
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 1]);
        }
      }
      if (iy == buffers->lens[0] - 1) {
        copy_center(&buffers->send_centers[6][ix], &centers[ind]);
        if(ix == buffers->lens[1] - 1) {
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 6]);
        }
      }
    }
  }
  for (long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    long base, stride;
    double x_offset, y_offset;
    switch (i) {
    case 0:
      base = 0;
      stride = 0;
      x_offset = (buffers->coords[1] == 0) ? -side : 0;
      y_offset = (buffers->coords[0] == 0) ? -side : 0;
      break;
    case 1:
      base = 1;
      stride = 1;
      x_offset = 0;
      y_offset = (buffers->coords[0] == 0) ? -side : 0;
      break;
    case 2:
      base = buffers->centers_lens[0] + 1;
      stride = 0;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = (buffers->coords[0] == 0) ? -side : 0;
      break;
    case 3:
      base = buffers->centers_lens[0] + 2;
      stride = buffers->centers_lens[0] + 2;
      x_offset = (buffers->coords[1] == 0) ? -side : 0;
      y_offset = 0;
      break;
    case 4:
      base = 2*buffers->centers_lens[0] + 1;
      stride = buffers->centers_lens[0] + 2;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = 0;
      break;
    case 5:
      base = (buffers->centers_lens[0] + 1)*(buffers->centers_lens[0] + 2);
      stride = 0;
      x_offset = (buffers->coords[1] == 0) ? -side : 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
      break;
    case 6:
      base = (buffers->centers_lens[0] + 1)*(buffers->centers_lens[0] + 2) + 1;
      stride = 1;
      x_offset = 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
      break;
    case 7:
      base = (buffers->centers_lens[0] + 2)*(buffers->centers_lens[0] + 2) - 1;
      stride = 0;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
      break;
    default:
      break;
    }
    MPI_Wait(&buffers->centers_requests[i], MPI_STATUS_IGNORE);
    for(long j = 0; j < buffers->centers_lens[i]; j++) {
      long ind = base + j*stride;
      centers[ind].x = buffers->recv_centers[i][j].x + x_offset;
      centers[ind].y = buffers->recv_centers[i][j].y + y_offset;
      centers[ind].m = buffers->recv_centers[i][j].m;
    }
    MPI_Start(&buffers->centers_requests[i]);
    MPI_Wait(&buffers->centers_requests[NUM_OF_NEIGHBORS + i], MPI_STATUS_IGNORE);
  }
}

void compute_kinetics(double side, long ncside, long ncells, cell_t *cells,
                      center_t *centers) {

  for (long i = 0; i < ncells; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      __m256d resx_vec = _mm256_setzero_pd();
      __m256d resy_vec = _mm256_setzero_pd();
      double mg = G * cells[i].m[j];
      __m256d mg2 = _mm256_set1_pd(mg);
      __m256d xj = _mm256_set1_pd(cells[i].x[j]);
      __m256d yj = _mm256_set1_pd(cells[i].y[j]);
      long long iter = (cell_size) - ((cell_size - (j + 1)) & 3);
      long long k;
      for (k = j + 1; k < iter; k += 4) {
        __m256d xk = _mm256_loadu_pd(&cells[i].x[k]);
        __m256d yk = _mm256_loadu_pd(&cells[i].y[k]);
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
        __m256d axk = _mm256_loadu_pd(&cells[i].ax[k]);
        __m256d ayk = _mm256_loadu_pd(&cells[i].ay[k]);
        axk = _mm256_sub_pd(axk, forcex);
        ayk = _mm256_sub_pd(ayk, forcey);
        _mm256_storeu_pd(&cells[i].ax[k], axk);
        _mm256_storeu_pd(&cells[i].ay[k], ayk);
        resx_vec = _mm256_add_pd(resx_vec, forcex);
        resy_vec = _mm256_add_pd(resy_vec, forcey);
      }
      double resx = sum_lanes(resx_vec);
      double resy = sum_lanes(resy_vec);
      while (k < cell_size) {
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

void compute_new_particle_cell(double size, long ncside, int id,
                               cell_t *cells,
                               communication_buffers_t *buffers) {
  for (long i = 0; i < buffers->size; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      int owner = find_owner_c(buffers, &cells[i], j, size, ncside);
      if (owner == id) {
        long ind = find_cell_c(buffers, &cells[i], j, size);
        cell_push_back_c(&cells[ind], &cells[i], j);
      } else {
        particles_buffer_add_c(&buffers->send_particles[owner], &cells[i], j);
      }
      cells[i].ind[j] = -1;
    }
  }
  for(long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    int neighbor = find_neighbor(buffers, i);
    MPI_Isend(buffers->send_particles[i].particles,
      buffers->send_particles[i].size * sizeof(particle_t), MPI_BYTE,
      neighbor, TAG_PARTICLE + NUM_OF_NEIGHBORS - i - 1, buffers->cart_comm,
      &buffers->particles_requests[i]);
    particles_buffer_resize(&buffers->send_particles[i], 0);
    while (!buffers->particles_flags[i]) {
      MPI_Iprobe(neighbor, TAG_PARTICLE + i, buffers->cart_comm,
                 &buffers->particles_flags[i], &buffers->particles_status[i]);
    }
    buffers->particles_flags[i] = 0;
    int bytes;
    MPI_Get_count(&buffers->particles_status[i], MPI_BYTE, &bytes);
    particles_buffer_resize(&buffers->recv_particles[i], bytes / sizeof(particle_t));
    MPI_Recv(buffers->recv_particles[i].particles, bytes, MPI_BYTE,
            neighbor, TAG_PARTICLE + i, buffers->cart_comm,
            &buffers->particles_status[i]);
    for (long long j = 0; j < buffers->recv_particles[i].size; j++) {
      long ind = find_cell_p(buffers, &buffers->recv_particles[i].particles[j], size);
      cell_push_back_p(&cells[ind], &buffers->recv_particles[i].particles[j]);
    }
    MPI_Iprobe(neighbor, TAG_PARTICLE + i, buffers->cart_comm,
             &buffers->particles_flags[i], &buffers->particles_status[i]);
    MPI_Wait(&buffers->particles_requests[i], MPI_STATUS_IGNORE);
  }
}

void detect_collisions(long ncells, cell_t *cells) {

  for (long i = 0; i < ncells; i++) {
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
      __m256d epsilon2 = _mm256_set1_pd(EPSILON2);
      __m256d xj = _mm256_set1_pd(cells[i].x[j]);
      __m256d yj = _mm256_set1_pd(cells[i].y[j]);
      long long k;
      long long iter = (j) - (j & 3);
      for (k = 0; k < iter; k += 4) {
        __m256d xk = _mm256_loadu_pd(&cells[i].x[k]);
        __m256d yk = _mm256_loadu_pd(&cells[i].y[k]);
        __m256d dx = _mm256_sub_pd(xk, xj);
        __m256d dy = _mm256_sub_pd(yk, yj);
        __m256d distance = _mm256_mul_pd(dx, dx);
        distance = _mm256_fmadd_pd(dy, dy, distance);
        __m256d cmp = _mm256_cmp_pd(distance, epsilon2, _CMP_LT_OQ);
        if (_mm256_testz_pd(cmp, cmp)) {
          continue;
        }
        __m256i near = _mm256_castpd_si256(cmp);
        __m256i collided = _mm256_loadu_si256((__m256i *)&cells[i].collided[k]);
        __m256i mask = _mm256_andnot_si256(collided, near);
        collision += (_mm256_movemask_epi8(mask) != 0);
        _mm256_storeu_si256((__m256i *)&cells[i].collided[k],
                            _mm256_or_si256(near, collided));
      }
      while (k < j) {
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
      if (collision > 0 && (cells[i].collided[j] == 0)) {
        cells[i].collided[j] = -1;
        cell_collisions++;
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
}

particle_t find_particle_zero(long id, long size, cell_t *cells) {
  particle_t par0;
  for (long i = 0; i < size; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      if (cells[i].ind[j] == 0) {
        par0.x = cells[i].x[j];
        par0.y = cells[i].y[j];
        par0.vx = cells[i].vx[j];
        par0.vy = cells[i].vy[j];
        par0.m = cells[i].m[j];
        par0.ind = 0;
        if(id != 0) {
          MPI_Send(&par0, sizeof(particle_t), MPI_BYTE, 0, TAG_PARTICLE_ZERO, MPI_COMM_WORLD);
        }
        return par0;
      }
    }
  }
  par0.x = -1;
  par0.y = -1;
  par0.vx = 0.0;
  par0.vy = 0.0;
  par0.m = 0.0;
  par0.ind = -1;
  return par0;
}

simulation_result simulation(double side, long ncside, long long npart, int id,
                             long long nstep, particles_buffer_t *par,
                             communication_buffers_t* buffers) {
  double size = side / ncside;
  cell_t *cells = malloc(sizeof(cell_t) * buffers->size);
  center_t *centers = malloc(sizeof(center_t) * ((buffers->lens[0] + 2)*(buffers->lens[1] + 2)));
  simulation_result res;
  init_block(size, ncside, npart, id, par, cells, buffers);
  for (long long i = 0; i < nstep; i++) {
    DEBUG("--------STEP %d: %lld --------------\n", id, i);
    compute_centers_of_mass(side, id, cells, centers, buffers);
    debug_centers(id, (buffers->lens[0] + 2)*(buffers->lens[1] + 2), centers);
    compute_kinetics(side, ncside, buffers->size, cells, centers);
    compute_new_particle_cell(size, ncside, id, cells, buffers);
    detect_collisions(buffers->size, cells);
    debug_particles(id, buffers->size, cells);
  }
  MPI_Reduce(&ncollisions, &res.number_of_collisions, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  res.particle_zero = find_particle_zero(id, buffers->size, cells);
  if(id == 0 && res.particle_zero.ind == -1) {
    MPI_Recv(&res.particle_zero, sizeof(particle_t), MPI_BYTE, MPI_ANY_SOURCE, TAG_PARTICLE_ZERO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  clean_blocks(cells, centers, buffers);
  return res;
}
