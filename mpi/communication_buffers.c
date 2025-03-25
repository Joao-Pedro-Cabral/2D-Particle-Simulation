
#include "communication_buffers.h"
#include "nodes.h"
#include <stdlib.h>

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
                                int id, int p, long long npart) {
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
    MPI_Recv_init(&buffers->recv_centers[i], sizeof(center_t) * len, MPI_BYTE, BLOCK_INDEX(id, i, p),
            TAG_CENTER + i, buffers->cart_comm, &buffers->centers_requests[i]);
    MPI_Send_init(&buffers->send_centers[i], sizeof(center_t) * len, MPI_BYTE, BLOCK_INDEX(id, i, p),
            TAG_CENTER + NUM_OF_NEIGHBORS - i - 1, buffers->cart_comm, &buffers->centers_requests[NUM_OF_NEIGHBORS + i]);
    MPI_Start(&buffers->centers_requests[i]);
    MPI_Iprobe(BLOCK_INDEX(id, i, p), TAG_PARTICLE + i, buffers->cart_comm,
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
}

int find_owner_xy(communication_buffers_t *buffers, double x, double y, double size, double ncside) {
  long xpart = x / size;
  long ypart = y / size;
  int row = xpart / (ncside / buffers->lens[0]);
  int col = ypart / (ncside / buffers->lens[1]);
  int coords[2] = {row, col};
  int owner;
  MPI_Cart_rank(buffers->cart_comm, coords, &owner);
  return owner;
}

int find_owner_c(communication_buffers_t *buffers, cell_t *cell, long long i, double size, double ncside) {
  long xpart = cell->x[i] / size;
  long ypart = cell->y[i] / size;
  int row = xpart / (ncside / buffers->lens[0]);
  int col = ypart / (ncside / buffers->lens[1]);
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
