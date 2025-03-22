#include "simulation.h"
#include "cells.h"
#include "debug.h"
#include "init_particles.h"
#include "particles.h"
#include "nodes.h"
#include <math.h>
#include <immintrin.h>
#include <mpi.h>
#include "communication_buffers.h"

static long long ncollisions = 0;

void init_blocks(double size, long ncside, long ncside2, long long npart,
                     int id, int p, particle_t *par, cell_t *cells,
                     communication_buffers_t *buffers) {
  long long *count = malloc(sizeof(long long) * BLOCK_SIZE(id, p, ncside));
  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
    count[i] = 0;
  }
  long block_low = BLOCK_LOW(id, p, ncside);
  long block_high = BLOCK_HIGH(id, p, ncside);
  for (long long i = 0; i < npart; i++) {
    long cell = find_cell_p(&par[i], size, ncside);
    if(block_low <= cell && cell <= block_high)
      count[cell - block_low]++;
  }
  long long min_size = 2*npart / ncside2;
  min_size = (min_size > 10) ? min_size : 10;
  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
    count[i] *= 2;
    long long capacity =
        count[i] > npart ? npart : (count[i] > min_size ? count[i] : min_size);
    cell_init(&cells[i], capacity,
              (i % ncside + 1) + (i / ncside + 1) * (ncside + 2));
  }
  for (long long i = 0; i < npart; i++) {
    long cell = find_cell_p(&par[i], size, ncside);
    if(block_low <= cell && cell <= block_high)
      cell_push_back_p(&cells[cell - block_low], &par[i]);
  }
  free(count);
  free(par);
  communication_buffers_init(buffers);
  MPI_Irecv(buffers->centers_up, sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id-1)%p, TAG_CENTER_UP, MPI_COMM_WORLD, &buffers->centers_requests[0]);
  MPI_Irecv(buffers->centers_down, sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id+1)%p, TAG_CENTER_DOWN, MPI_COMM_WORLD, &buffers->centers_requests[1]);
  MPI_Iprobe((id-1)%p, TAG_PARTICLE_UP, MPI_COMM_WORLD, &buffers->particles_flags[0], &buffers->particles_status[0]);
  MPI_Iprobe((id+1)%p, TAG_PARTICLE_DOWN, MPI_COMM_WORLD, &buffers->particles_flags[1], &buffers->particles_status[1]);
}

void clean_blocks(long ncside, int id, int p, cell_t *cells,
                  center_t *centers, communication_buffers_t *buffers) {
  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
    cell_clean(&cells[i]);
  }
  free(cells);
  free(centers);
  communication_buffers_clean(&buffers);
}

void debug_centers(long block_size, center_t *centers) {
  for (long i = 0; i < block_size; i++) {
    DEBUG("Center %ld x: %.6lf y: %.6lf m: %.6lf\n", i, centers[i].x,
          centers[i].y, centers[i].m);
  }
}

static inline double sum_lanes(__m256d vec) {
  vec = _mm256_hadd_pd(vec, vec);
  return ((double*)&vec)[0] + ((double*)&vec)[2];
}

void compute_centers_of_mass(double side, long ncside, int id, int p,
                             cell_t *cells, center_t *centers,
                             communication_buffers_t * buffers) {

  for (long iy = 0; iy < NUM_ROWS(id, p, ncside); iy++) {
    for (long ix = 0; ix < NUM_COLUMNS(id, p, ncside); ix++) {
    long i = ix + ncside*iy;
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

    if(ix == 0 || ix == NUM_COLUMNS(id, p, ncside) - 1) {
      long ind2 = (ix == 0) ? ind + NUM_COLUMNS(id, p, ncside) : ind - NUM_COLUMNS(id, p, ncside);
      centers[ind2].x = (ix == 0) ? centers[ind].x + side : centers[ind].x - side;
      centers[ind2].y = centers[ind].y;
      centers[ind2].m = centers[ind].m;
    }
  }
  }
  MPI_Isend(&centers[cells[0].center], sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id-1)%p, TAG_CENTER_DOWN, MPI_COMM_WORLD, &buffers->centers_requests[2]);
  MPI_Isend(&centers[cells[ncside*(NUM_ROWS(id,p,ncside)-1)].center], sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id+1)%p, TAG_CENTER_UP, MPI_COMM_WORLD, &buffers->centers_requests[3]);
  MPI_Waitall(2*BLOCK_NUM_OF_NEIGHBORS(id,p,ncside), buffers->centers_requests, MPI_STATUS_IGNORE);
  for (long i = 0; i < NUM_COLUMNS(id, p, ncside); i++) {
    long ind = i + 1;
    centers[ind].x = buffers->centers_up[i].x;
    centers[ind].y = (id == 0) ? buffers->centers_up[i].y - side : buffers->centers_up[i].y;
    centers[ind].m = buffers->centers_up[i].m;
    long ind2 = NUM_ROWS(id, p, ncside)*NUM_COLUMNS(id, p, ncside) + ind;
    centers[ind2].x = buffers->centers_down[i].x;
    centers[ind2].y = (id == p - 1) ? buffers->centers_down[i].y + side : buffers->centers_down[i].y;
    centers[ind2].m = buffers->centers_down[i].m;
    if(i == 0 || i == NUM_COLUMNS(id, p, ncside) - 1) {
      long ind2 = (i == 0) ? ind + NUM_COLUMNS(id, p, ncside) : ind - NUM_COLUMNS(id, p, ncside);
      centers[ind2].x = (i == 0) ? centers[ind].x + side : centers[ind].x - side;
      centers[ind2].y = centers[ind].y;
      centers[ind2].m = centers[ind].m;
    }
  }
  MPI_Irecv(buffers->centers_up, sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id-1)%p, TAG_CENTER_UP, MPI_COMM_WORLD, &buffers->centers_requests[0]);
  MPI_Irecv(buffers->centers_down, sizeof(center_t)*NUM_COLUMNS(id, p, ncside), MPI_BYTE, (id+1)%p, TAG_CENTER_DOWN, MPI_COMM_WORLD, &buffers->centers_requests[1]);
  debug_centers(BLOCK_SIZE(id, p, ncside) + BLOCK_NEIGHBORHOOD(id, p, ncside), centers);
}

void compute_kinetics(double side, long ncside, int id, int p, cell_t *cells,
                      center_t *centers) {

  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
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

void debug_particles(long block_size, cell_t *cells) {
  for (long i = 0; i < block_size; i++) {
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

void compute_new_particle_cell(double size, long ncside, int id, int p, cell_t *cells, 
                               communication_buffers_t * buffers) {
  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      long ind = find_cell_c(&cells[i], j, size, ncside);
      if (ind == i)
        continue;
      int owner = BLOCK_OWNER(ind, p, n);
      if(owner == id) {
        cell_push_back_c(&cells[ind], &cells[i], j);
      } else if (owner == (id - 1) % p) {
        particles_buffer_add(&buffers->send_particles_up, &cells[i], j);
      } else if (owner == (id + 1) % p) {
        particles_buffer_add(&buffers->send_particles_down, &cells[i], j);
      }
      cells[i].ind[j] = -1;
    }
  }
  MPI_Send(&buffers->send_particles_up, buffers->send_particles_up.size*sizeof(particle_t), MPI_BYTE, (id - 1) % p, TAG_PARTICLE_UP, MPI_COMM_WORLD);
  MPI_Send(&buffers->send_particles_down, buffers->send_particles_down.size*sizeof(particle_t), MPI_BYTE, (id + 1) % p, TAG_PARTICLE_DOWN, MPI_COMM_WORLD);
  while (!buffers->particles_flags[0] || !buffers->particles_flags[1]) {
    MPI_Iprobe((id-1)%p, TAG_PARTICLE_UP, MPI_COMM_WORLD, &buffers->particles_flags[0], &buffers->particles_status[0]);
    MPI_Iprobe((id+1)%p, TAG_PARTICLE_DOWN, MPI_COMM_WORLD, &buffers->particles_flags[1], &buffers->particles_status[1]);
  }
  int bytes_up, bytes_down;
  MPI_Get_count(&buffers->particles_status[0], MPI_BYTE, &bytes_up);
  MPI_Get_count(&buffers->particles_status[1], MPI_BYTE, &bytes_down);
  particles_buffer_resize(&buffers->recv_particles_up, bytes_up/sizeof(particle_t));
  particles_buffer_resize(&buffers->recv_particles_down, bytes_down/sizeof(particle_t));
  MPI_Recv(&buffers->recv_particles_up, bytes_up, MPI_BYTE, (id - 1) % p, TAG_PARTICLE_UP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&buffers->recv_particles_down, bytes_down, MPI_BYTE, (id + 1) % p, TAG_PARTICLE_DOWN,, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  for (long i = 0; i < bytes_up/sizeof(particle_t); i++) {
    long ind = find_cell_p(&buffers->recv_particles_up[i], size, ncside);
    cell_push_back_p(&cells[ind], &buffers->recv_particles_up[i]);
  }
  for (long i = 0; i < bytes_down/sizeof(particle_t); i++) {
    long ind = find_cell_p(&buffers->recv_particles_down[i], size, ncside);
    cell_push_back_p(&cells[ind], &buffers->recv_particles_down[i]);
  }
}

void detect_collisions(long ncside, int id, int p, cell_t *cells) {

  for (long i = 0; i < BLOCK_SIZE(id, p, ncside); i++) {
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
  par0.x = -1;
  par0.y = -1;
  par0.vx = 0.0;
  par0.vy = 0.0;
  par0.m = 0.0;
  return par0;
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep, particle_t *par) {
  double size = side / ncside;
  long ncside2 = ncside * ncside;
  int id, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  cell_t *cells = malloc(sizeof(cell_t) * BLOCK_SIZE(id, p, ncside));
  center_t *centers = malloc(sizeof(center_t) * (BLOCK_SIZE(id, p, ncside) + BLOCK_NEIGHBORHOOD(id, p, ncside)));
  simulation_result res;
  communication_buffers_t buffers;
  init_blocks(size, ncside, ncside2, npart, id, p, par, cells, &buffers);
  for (long long i = 0; i < nstep; i++) {
    // printf("--------STEP: %lld --------------\n", i);
    compute_centers_of_mass(side, ncside, id, p, cells, centers, &buffers);
    compute_kinetics(side, ncside, id, p, cells, centers);
    compute_new_particle_cell(size, ncside,  id, p, cells, &buffers);
    detect_collisions(ncside, id, p, cells);
  }
  res.number_of_collisions = ncollisions;
  res.particle_zero = find_particle_zero(ncside2, cells);
  clean_blocks(ncside, id, p, cells, centers, &buffers);
  return res;
}
