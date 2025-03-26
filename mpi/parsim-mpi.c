#include "debug.h"
#include "init_particles.h"
#include "nodes.h"
#include "simulation.h"
#include "utils.h"
#include "particles_buffer.h"
#include "communication_buffers.h"
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void print_result(simulation_result res) {
  printf("%.3lf %.3lf\n", res.particle_zero.x, res.particle_zero.y);
  printf("%lld\n", res.number_of_collisions);
}

int main(int argc, char *argv[]) {
  DEBUG("%d\n", argc);
  if (argc != 6) {
    ERROR("Incorrect number of arguments\n");
  }

  long seed;
  double side;
  long ncside;
  long long npart, nstep;

  if (!cast_to_long(argv[1], &seed)) {
    ERROR("Unexpected format for seed\n");
  }

  if (!cast_to_double(argv[2], &side)) {
    ERROR("Unexpected format for side\n");
  }

  if (!cast_to_long(argv[3], &ncside)) {
    ERROR("Unexpected format for ncside\n");
  }

  if (!cast_to_long_long(argv[4], &npart)) {
    ERROR("Unexpected format for npart\n");
  }

  if (!cast_to_long_long(argv[5], &nstep)) {
    ERROR("Unexpected format for nstep\n");
  }

  DEBUG("%ld, %lf, %ld, %lld, %lld\n", seed, side, ncside, npart, nstep);

  MPI_Init(&argc, &argv);

  int id, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  communication_buffers_t buffers;
  communication_buffers_create_world(&buffers, ncside, id, p);

  particles_buffer_t par;
  particles_buffer_init(&par, npart/p);

  double exec_time;
  init_particles(seed, side, ncside, id, npart, &par, &buffers);
  exec_time = -omp_get_wtime();
  simulation_result res = simulation(side, ncside, npart, id, nstep, &par, &buffers);
  exec_time += omp_get_wtime();
  if(id == 0) {
    print_result(res);
    fprintf(stderr, "%.1fs\n", exec_time);
  }
  fflush(stdout);
  MPI_Finalize();
  return 0;
}
