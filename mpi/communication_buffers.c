
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
    int neighbor = find_neighbor(buffers, id, i);
    DEBUG("Process: %d, Neighbor: %d, len: %ld, recv: %ld, send: %ld\n", id, neighbor, len, TAG_CENTER + i, TAG_CENTER + NUM_OF_NEIGHBORS - i - 1);
    MPI_Recv_init(&buffers->recv_centers[i], sizeof(center_t) * len, MPI_BYTE, neighbor,
            TAG_CENTER + i, buffers->cart_comm, &buffers->centers_requests[i]);
    MPI_Send_init(&buffers->send_centers[i], sizeof(center_t) * len, MPI_BYTE, neighbor,
            TAG_CENTER + NUM_OF_NEIGHBORS - i - 1, buffers->cart_comm, &buffers->centers_requests[NUM_OF_NEIGHBORS + i]);
    MPI_Start(&buffers->centers_requests[i]);
    MPI_Iprobe(neighbor, TAG_PARTICLE + i, buffers->cart_comm,
             &buffers->particles_flags[i], &buffers->particles_status[i]); // TODO: Good idea?
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

int find_owner_xy(communication_buffers_t *buffers, double x, double y, double size, long ncside) {
  long xpart = x / size;
  long ypart = y / size;
  int col = (buffers->dims[1]*(xpart+1)-1) / ncside;
  int row = (buffers->dims[0]*(ypart+1)-1) / ncside;
  int coords[2] = {row, col};
  int owner;
  MPI_Cart_rank(buffers->cart_comm, coords, &owner);
  return owner;
}

int find_owner_c(communication_buffers_t *buffers, cell_t *cell, long long i, double size, long ncside) {
  long xpart = cell->x[i] / size;
  long ypart = cell->y[i] / size;
  int col = (buffers->dims[1]*(xpart+1)-1) / ncside;
  int row = (buffers->dims[0]*(ypart+1)-1) / ncside;
  int coords[2] = {row, col};
  int owner;
  MPI_Cart_rank(buffers->cart_comm, coords, &owner);
  return owner;
}

long find_cell_p(communication_buffers_t *buffers, particle_t *par, double size) {
  long xpart = par->x / size;
  long ypart = par->y / size;
  return buffers->lens[0] * (ypart - buffers->rows[0]) + (xpart - buffers->cols[0]);
}

long find_cell_c(communication_buffers_t *buffers, cell_t *cell, long long i, double size) {
  long xpart = cell->x[i] / size;
  long ypart = cell->y[i] / size;
  return buffers->lens[0] * (ypart - buffers->rows[0]) + (xpart - buffers->cols[0]);
}

int find_neighbor(communication_buffers_t *buffers, int id, int pos) {
  int neighbor, aux;
  switch (pos) {
  case 0:
    MPI_Cart_shift(buffers->cart_comm, 1, -1, &id, &aux);
    MPI_Cart_shift(buffers->cart_comm, 0, -1, &aux, &neighbor);
    break;
  case 1:
    MPI_Cart_shift(buffers->cart_comm, 0, -1, &id, &neighbor);
    break;
  case 2:
    MPI_Cart_shift(buffers->cart_comm, 0, -1, &id, &aux);
    MPI_Cart_shift(buffers->cart_comm, 1, 1, &aux, &neighbor);
    break;
  case 3:
    MPI_Cart_shift(buffers->cart_comm, 1, -1, &id, &neighbor);
    break;
  case 4:
    MPI_Cart_shift(buffers->cart_comm, 1, 1, &id, &neighbor);
    break;
  case 5:
    MPI_Cart_shift(buffers->cart_comm, 0, 1, &id, &aux);
    MPI_Cart_shift(buffers->cart_comm, 1, -1, &aux, &neighbor);
    break;
  case 6:
    MPI_Cart_shift(buffers->cart_comm, 0, 1, &id, &neighbor);
    break;
  default: // 7
    MPI_Cart_shift(buffers->cart_comm, 1, 1, &id, &aux);
    MPI_Cart_shift(buffers->cart_comm, 0, 1, &aux, &neighbor);
    break;
  }
  return neighbor;
}
