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
    DEBUG("From %d: Center %ld x: %.15lf y: %.15lf m: %.15lf\n", id, i,
          centers[i].x, centers[i].y, centers[i].m);
  }
}

void debug_particles(int id, long ncells, cell_t *cells) {
  for (long i = 0; i < ncells; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      DEBUG("From %d: Particle %lld: m: %.15f, x: %.15f, y: %.15f, vx: %.15f, vy: "
            "%.15f, "
            "ax: %.15f, ay: %.15f, collided: %lld, cell: %ld\n",
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

void compute_centers_of_mass(double side,
                             cell_t *cells, center_t *centers,
                             communication_buffers_t *buffers,
                             long chunk_size) {

  #pragma omp for collapse(2) schedule(dynamic, chunk_size)
  for (long iy = 0; iy < buffers->lens[0]; iy++) {
    for (long ix = 0; ix < buffers->lens[1]; ix++) {
      long i = ix + buffers->lens[1] * iy;
      double total_mass = 0.0;
      double weighted_x = 0.0;
      double weighted_y = 0.0;
      long long cell_size = cells[i].size;

      long long j;
      for (j = 0; j < cell_size; j ++) {
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

      if (ix == 0) {
        copy_center(&buffers->send_centers[3][iy], &centers[ind]);
        if(iy == 0) {
          copy_center(&buffers->send_centers[0][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 0]);
        }
        if (iy == buffers->lens[0] - 1) {
          copy_center(&buffers->send_centers[5][0], &centers[ind]);
          MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 5]);
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
        }
      }
      if(iy == 0) {
        copy_center(&buffers->send_centers[1][ix], &centers[ind]);
      }
      if (iy == buffers->lens[0] - 1) {
        copy_center(&buffers->send_centers[6][ix], &centers[ind]);
      }
    }
  }
  #pragma omp single nowait
  {
  MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 1]);
  MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 3]);
  MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 4]);
  MPI_Start(&buffers->centers_requests[NUM_OF_NEIGHBORS + 6]);
  }
  #pragma omp for schedule(dynamic, 1)
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
      base = buffers->lens[1] + 1;
      stride = 0;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = (buffers->coords[0] == 0) ? -side : 0;
      break;
    case 3:
      base = buffers->lens[1] + 2;
      stride = buffers->lens[1] + 2;
      x_offset = (buffers->coords[1] == 0) ? -side : 0;
      y_offset = 0;
      break;
    case 4:
      base = 2*buffers->lens[1] + 3;
      stride = buffers->lens[1] + 2;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = 0;
      break;
    case 5:
      base = (buffers->lens[0] + 1)*(buffers->lens[1] + 2);
      stride = 0;
      x_offset = (buffers->coords[1] == 0) ? -side : 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
      break;
    case 6:
      base = (buffers->lens[0] + 1)*(buffers->lens[1] + 2) + 1;
      stride = 1;
      x_offset = 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
      break;
    default: // 7
      base = (buffers->lens[0] + 2)*(buffers->lens[1] + 2) - 1;
      stride = 0;
      x_offset = (buffers->coords[1] == buffers->dims[1] - 1) ? side : 0;
      y_offset = (buffers->coords[0] == buffers->dims[0] - 1) ? side : 0;
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
  }
}


void compute_kinetics(double side, cell_t *cells,
                      center_t *centers, const communication_buffers_t *buffers,
                      long chunk_size) {

  #pragma omp for schedule(dynamic, chunk_size)
  for (long i = 0; i < buffers->size; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      double mg = G * cells[i].m[j];
      double resx = 0.0;
      double resy = 0.0;     
      long long k;
      for (k = j + 1; k < cell_size; k ++) {
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
        long ind = cells[i].center + ((k % 3) - 1) + (k / 3 - 1) * (buffers->lens[1] + 2);
        gravitational_force_pc(&cells[i], j, &centers[ind]);
      }
      cells[i].ax[j] /= cells[i].m[j];
      cells[i].ay[j] /= cells[i].m[j];
      update_position_and_velocity(&cells[i], j, side);
    }
  }
}

void compute_new_particle_cell(double size, long ncside, int id,
                               cell_t *cells, communication_buffers_t *buffers,
                               long chunk_size) {
  for (long i = 0; i < buffers->size; i++) {
    long long cell_size = cells[i].size;
    for (long long j = 0; j < cell_size; j++) {
      int owner = find_owner_c(buffers, &cells[i], j, size, ncside);
      if (owner == id) {
        long ind = find_cell_c(buffers, &cells[i], j, size);
        if(ind == i)
          continue;
        cell_push_back_c(&cells[ind], &cells[i], j);
      } else {
        int pos = find_neighbor_pos(buffers, owner);
        particles_buffer_add_c(&buffers->send_particles[pos], &cells[i], j);
      }
      cells[i].ind[j] = -1;
    }
  }
  for(long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    int neighbor = buffers->neighbors[i];
    MPI_Isend(buffers->send_particles[i].particles,
      buffers->send_particles[i].size * sizeof(particle_t), MPI_BYTE,
      neighbor, TAG_PARTICLE + NUM_OF_NEIGHBORS - i - 1, buffers->cart_comm,
      &buffers->particles_requests[i]);
    particles_buffer_resize(&buffers->send_particles[i], 0);
  }
  for(long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    int neighbor = buffers->neighbors[i];
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
  }
  MPI_Waitall(NUM_OF_NEIGHBORS, &buffers->centers_requests[NUM_OF_NEIGHBORS], MPI_STATUSES_IGNORE);
  MPI_Waitall(NUM_OF_NEIGHBORS, &buffers->particles_requests[0], MPI_STATUSES_IGNORE);
}

void detect_collisions(long ncells, cell_t *cells, long chunk_size) {

  #pragma omp for reduction(+ : ncollisions) schedule(dynamic, chunk_size)
  for (long i = 0; i < ncells; i++) {
    long long cell_collisions = 0;
    long long cell_size = cells[i].size;
    for (long long j = cell_size - 1; j >= 0; j--) {
      if (cells[i].ind[j] == -1) {
        // DEBUG("Removing mass null Particle %lld from cells[%ld].par[%lld]\n",
        //       cells[i].ind[j], i, j);
        cell_size--;
        cell_remove_particle(&cells[i], j);
      }
    }
    for (long long j = cell_size - 1; j >= 0; j--) {
      long long collision = 0;
      long long k;
      for (k = 0; k < j; k ++) {
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
  long chunk_size = (buffers->size >= 4*omp_get_max_threads()) ? 2 : 1;
  #pragma omp parallel
  for (long long i = 0; i < nstep; i++) {
    // DEBUG("--------STEP %d: %lld --------------\n", id, i);
    compute_centers_of_mass(side, cells, centers, buffers, chunk_size);
    compute_kinetics(side, cells, centers, buffers, chunk_size);
    #pragma omp single
    {
    compute_new_particle_cell(size, ncside, id, cells, buffers, chunk_size);
    }
    detect_collisions(buffers->size, cells, chunk_size);
  }
  MPI_Reduce(&ncollisions, &res.number_of_collisions, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  res.particle_zero = find_particle_zero(id, buffers->size, cells);
  if(id == 0 && res.particle_zero.ind == -1) {
    MPI_Recv(&res.particle_zero, sizeof(particle_t), MPI_BYTE, MPI_ANY_SOURCE, TAG_PARTICLE_ZERO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  clean_blocks(cells, centers, buffers);
  return res;
}