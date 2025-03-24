
#include "communication_buffers.h"
#include "nodes.h"
#include <stdlib.h>

void communication_buffers_init(communication_buffers_t *buffers, long ncside,
                                int id, int p, long long npart) {
  buffers->send_centers_up =
      malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->send_centers_down =
      malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->recv_centers_up =
      malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->recv_centers_down =
      malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->centers_requests =
      malloc(sizeof(MPI_Request) * 2 * BLOCK_NUM_OF_NEIGHBORS(id, p, ncside));
  long long initial_estimation = 2 * npart / ncside;
  particles_buffer_init(&buffers->recv_particles_up, initial_estimation);
  particles_buffer_init(&buffers->recv_particles_down, initial_estimation);
  particles_buffer_init(&buffers->send_particles_up, initial_estimation);
  particles_buffer_init(&buffers->send_particles_down, initial_estimation);
  buffers->particles_flags =
      malloc(sizeof(int) * BLOCK_NUM_OF_NEIGHBORS(id, p, ncside));
  for (int i = 0; i < BLOCK_NUM_OF_NEIGHBORS(id, p, ncside); i++) {
    buffers->particles_flags[i] = 0;
  }
  buffers->particles_status =
      malloc(sizeof(MPI_Status) * BLOCK_NUM_OF_NEIGHBORS(id, p, ncside));
  buffers->particles_requests = malloc(sizeof(MPI_Request) * BLOCK_NUM_OF_NEIGHBORS(id, p, ncside));
}

void communication_buffers_clean(communication_buffers_t *buffers) {
  free(buffers->send_centers_up);
  free(buffers->send_centers_down);
  free(buffers->recv_centers_up);
  free(buffers->recv_centers_down);
  free(buffers->centers_requests);
  particles_buffer_clean(&buffers->recv_particles_up);
  particles_buffer_clean(&buffers->recv_particles_down);
  particles_buffer_clean(&buffers->send_particles_up);
  particles_buffer_clean(&buffers->send_particles_down);
  free(buffers->particles_status);
  free(buffers->particles_flags);
  free(buffers->particles_requests);
  buffers->send_centers_up = NULL;
  buffers->send_centers_down = NULL;
  buffers->recv_centers_up = NULL;
  buffers->recv_centers_down = NULL;
  buffers->centers_requests = NULL;
  buffers->particles_flags = NULL;
  buffers->particles_status = NULL;
  buffers->particles_requests = NULL;
}
