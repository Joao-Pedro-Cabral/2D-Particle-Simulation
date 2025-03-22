
#include "communication_buffers.h"
#include "nodes.h"

void communication_buffers_init(communication_buffers_t *buffers, long ncside, int id, int p, long long npart) {
  buffers->centers_up = malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->centers_down = malloc(sizeof(center_t) * (BLOCK_FRONTIER(id, p, ncside)));
  buffers->centers_requests = malloc(sizeof(MPI_Request)*2*BLOCK_NUM_OF_NEIGHBORS(id,p,ncside));
  particles_buffer_init(buffers->recv_particles_up);
  particles_buffer_init(buffers->recv_particles_down);
  particles_buffer_init(buffers->send_particles_up);
  particles_buffer_init(buffers->send_particles_down);
  buffers->particles_flags = malloc(sizeof(int)*BLOCK_NUM_OF_NEIGHBORS(id,p,ncside));
  buffers->particles_status = malloc(sizeof(MPI_Status)*BLOCK_NUM_OF_NEIGHBORS(id,p,ncside));
}

void communication_buffers_clean(communication_buffers_t *buffers) {
  free(buffers->centers_up);
  free(buffers->centers_down);
  free(buffers->centers_requests);
  particles_buffer_clean(buffers->recv_particles_up);
  particles_buffer_clean(buffers->recv_particles_down);
  particles_buffer_clean(buffers->send_particles_up);
  particles_buffer_clean(buffers->send_particles_down);
  free(buffers->particles_status);
  free(buffers->particles_flags);
  buffers->centers_up = NULL;
  buffers->centers_down = NULL;
  buffers->centers_requests = NULL;
  buffers->particles_flags = NULL;
  buffers->particles_status = NULL;
}
