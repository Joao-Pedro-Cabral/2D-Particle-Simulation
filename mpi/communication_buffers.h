
#ifndef __COMMUNICATION_BUFFERS_H
#define __COMMUNICATION_BUFFERS_H

#include "particles_buffer.h"
#include "nodes.h"
#include <mpi.h>

typedef struct {
  MPI_Comm cart_comm;
  int dims [2];
  int coords [2];
  long rows [2];
  long cols [2];
  long lens [2];
  long size;
  int centers_lens [NUM_OF_NEIGHBORS];
  center_t *send_centers [NUM_OF_NEIGHBORS];
  center_t *recv_centers [NUM_OF_NEIGHBORS];
  MPI_Request centers_requests [2*NUM_OF_NEIGHBORS];
  particles_buffer_t send_particles [NUM_OF_NEIGHBORS];
  particles_buffer_t recv_particles [NUM_OF_NEIGHBORS];
  int particles_flags [NUM_OF_NEIGHBORS];
  MPI_Status particles_status [NUM_OF_NEIGHBORS];
  MPI_Request particles_requests [NUM_OF_NEIGHBORS];
} communication_buffers_t;

void communication_buffers_create_world(communication_buffers_t *buffers, long ncside, int id, int p);
void communication_buffers_init(communication_buffers_t *buffers, long ncside, int id, long long npart);
void communication_buffers_clean(communication_buffers_t *buffers);
int find_owner_xy(communication_buffers_t *buffers, double x, double y, double size, long ncside);
int find_owner_c(communication_buffers_t *buffers, cell_t * cell, long long i, double size, long ncside);
long find_cell_p(communication_buffers_t *buffers, particle_t *par, double size);
long find_cell_c(communication_buffers_t *buffers, cell_t * cell, long long i, double size);
int find_neighbor(communication_buffers_t *buffers, int pos);

#endif // __COMMUNICATION_BUFFERS_H