#include "debug.h"
#include "init_particles.h"
#include "parsim.h"
#include "utils.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

void print_result(const simulation_result &res) {
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

  std::vector<particle_t> par(npart);

  double exec_time;
  init_particles(seed, side, ncside, npart,
                 par.data()); // .data() to provide C compatibility
  exec_time = -omp_get_wtime();
  simulation_result res = simulation(side, ncside, npart, nstep, par);
  exec_time += omp_get_wtime();
  fprintf(stderr, "%.1fs\n", exec_time);
  print_result(res);

  return 0;
}
