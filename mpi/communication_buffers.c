
#include "communication_buffers.h"
#include "nodes.h"
#include <stdlib.h>
#include <stdio.h>
#include "debug.h"

void communication_buffers_create_world(communication_buffers_t *buffers, long ncside, int id, int p) {
  buffers->dims[0] = buffers->dims[1] = 0;
  MPI_Dims_create(p, 2, buffers->dims);
  int periods[2] = {1, 1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, buffers->dims, periods, 0, &buffers->cart_comm);
  MPI_Cart_coords(buffers->cart_comm, id, 2, buffers->coords);
  buffers->rows[0] = (buffers->coords[0] * ncside) / buffers->dims[0];
  buffers->rows[1] = ((buffers->coords[0] + 1) * ncside) / buffers->dims[0];
  buffers->cols[0] = (buffers->coords[1] * ncside) / buffers->dims[1];
  buffers->cols[1] = ((buffers->coords[1] + 1) * ncside) / buffers->dims[1];
  buffers->lens[0] = buffers->rows[1] - buffers->rows[0];
  buffers->lens[1] = buffers->cols[1] - buffers->cols[0];
  buffers->size = buffers->lens[0]*buffers->lens[1];
}

void communication_buffers_init(communication_buffers_t *buffers, long ncside,
                                int id, long long npart) {
  long long initial_estimation = 2 * npart / ncside;
  DEBUG("Process: %d, lens: %ld, %ld, row: %ld, col: %ld, dims: %d, %d\n", id, buffers->lens[0], buffers->lens[1], buffers->rows[0], buffers->cols[0], buffers->dims[0], buffers->dims[1]);
  for(long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    long len = 1;
    if (i == 1 || i == 6) {
      len = buffers->lens[1];
    } else if (i == 3 || i == 4) {
      len = buffers->lens[0];
    }
    buffers->centers_lens[i] = len;
    buffers->send_centers[i] = (center_t*) malloc(sizeof(center_t)*len);
    buffers->recv_centers[i] = (center_t*) malloc(sizeof(center_t)*len);
    particles_buffer_init(&buffers->send_particles[i], initial_estimation);
    particles_buffer_init(&buffers->recv_particles[i], initial_estimation);
    buffers->particles_flags[i] = 0;
    buffers->neighbors[i] = find_neighbor(buffers, i);
    DEBUG("Process: %d, Neighbor: %d, len: %ld, recv: %ld, send: %ld\n", id, buffers->neighbors[i], len, TAG_CENTER + i, TAG_CENTER + NUM_OF_NEIGHBORS - i - 1);
    MPI_Recv_init(buffers->recv_centers[i], sizeof(center_t) * len, MPI_BYTE, buffers->neighbors[i],
            TAG_CENTER + i, buffers->cart_comm, &buffers->centers_requests[i]);
    MPI_Send_init(buffers->send_centers[i], sizeof(center_t) * len, MPI_BYTE, buffers->neighbors[i],
            TAG_CENTER + NUM_OF_NEIGHBORS - i - 1, buffers->cart_comm, &buffers->centers_requests[NUM_OF_NEIGHBORS + i]);
    MPI_Start(&buffers->centers_requests[i]);
  }
}

void communication_buffers_clean(communication_buffers_t *buffers) {
  for(long i = 0; i < NUM_OF_NEIGHBORS; i++) {
    free(buffers->send_centers[i]);
    free(buffers->recv_centers[i]);
    buffers->send_centers[i] = NULL;
    buffers->recv_centers[i] = NULL;
    particles_buffer_clean(&buffers->recv_particles[i]);
    particles_buffer_clean(&buffers->send_particles[i]);
  }
  MPI_Comm_free(&buffers->cart_comm);
}

int find_owner_c(communication_buffers_t *buffers, cell_t *cell, long long i, double size, long ncside) {
  long xpart = cell->x[i] / size;
  long ypart = cell->y[i] / size;
  int col = (buffers->dims[1]*(xpart+1)-1) / ncside;
  int row = (buffers->dims[0]*(ypart+1)-1) / ncside;
  return col + row*buffers->dims[1];
}

long find_cell_p(communication_buffers_t *buffers, particle_t *par, double size) {
  long xpart = par->x / size;
  long ypart = par->y / size;
  return buffers->lens[1] * (ypart - buffers->rows[0]) + (xpart - buffers->cols[0]);
}

long find_cell_c(communication_buffers_t *buffers, cell_t *cell, long long i, double size) {
  long xpart = cell->x[i] / size;
  long ypart = cell->y[i] / size;
  return buffers->lens[1] * (ypart - buffers->rows[0]) + (xpart - buffers->cols[0]);
}

int find_neighbor(communication_buffers_t *buffers, int pos) {
  int neighbor, coords[2];
  switch (pos) {
  case 0:
    coords[0] = buffers->coords[0] - 1;
    coords[1] = buffers->coords[1] - 1;
    break;
  case 1:
    coords[0] = buffers->coords[0] - 1;
    coords[1] = buffers->coords[1];
    break;
  case 2:
    coords[0] = buffers->coords[0] - 1;
    coords[1] = buffers->coords[1] + 1;
    break;
  case 3:
    coords[0] = buffers->coords[0];
    coords[1] = buffers->coords[1] - 1;
    break;
  case 4:
    coords[0] = buffers->coords[0];
    coords[1] = buffers->coords[1] + 1;
    break;
  case 5:
    coords[0] = buffers->coords[0] + 1;
    coords[1] = buffers->coords[1] - 1;
    break;
  case 6:
    coords[0] = buffers->coords[0] + 1;
    coords[1] = buffers->coords[1];
    break;
  default: // 7
    coords[0] = buffers->coords[0] + 1;
    coords[1] = buffers->coords[1] + 1;
    break;
  }
  MPI_Cart_rank(buffers->cart_comm, coords, &neighbor);
  return neighbor;
}

int find_neighbor_pos(communication_buffers_t *buffers, int neighbor) {
  for(int i = 0; i < NUM_OF_NEIGHBORS; i++) {
    if(neighbor == buffers->neighbors[i]) {
      return i;
    }
  }
  return -1;
}
