
#ifndef __COMMUNICATION_BUFFERS_H
#define __COMMUNICATION_BUFFERS_H

#include "particles_buffer.h"
#include <mpi.h>

typedef struct {
  center_t *centers_up;
  center_t* centers_down;
  particles_buffer_t send_particles_up;
  particles_buffer_t send_particles_down;
  particles_buffer_t recv_particles_up;
  particles_buffer_t recv_particles_down;
  MPI_Request * centers_requests;
  int * particles_flags;
  MPI_Status * particles_status;
} communication_buffers_t;

void communication_buffers_init(communication_buffers_t *buffers, long ncside, int id, int p, long long npart);
void communication_buffers_clean(communication_buffers_t *buffers);

#endif // __COMMUNICATION_BUFFERS_H