#include <stdio.h>
#include "debug.h"
#include "utils.h"
#include "init_particles.h"
#include <omp.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
  DEBUG("%d\n", argc);
  if(argc != 6) {
    DEBUG("Incorrect number of arguments\n");
    return 1;
  }

  long seed;
  double side;
  long ncside;
  long long npart, nstep;

  if(!cast_to_long(argv[1], &seed)) {
    DEBUG("Unexpected format for seed\n");
    return 1;
  }

  if(!cast_to_double(argv[2], &side)) {
    DEBUG("Unexpected format for side\n");
    return 1;
  }

  if(!cast_to_long(argv[3], &ncside)) {
    DEBUG("Unexpected format for ncside\n");
    return 1;
  }

  if(!cast_to_long_long(argv[4], &npart)) {
    DEBUG("Unexpected format for npart\n");
    return 1;
  }

  if(!cast_to_long_long(argv[5], &nstep)) {
    DEBUG("Unexpected format for nstep\n");
    return 1;
  }

  DEBUG("%ld, ", seed);
  DEBUG("%lf, ", side);
  DEBUG("%ld, ", ncside);
  DEBUG("%lld, ", npart);
  DEBUG("%lld\n", nstep);

  particle_t* par = malloc(npart*sizeof(particle_t));


  double exec_time;
  init_particles(seed, side, ncside, npart, par);
  exec_time =-omp_get_wtime();
  //simulation();
  exec_time += omp_get_wtime();
  fprintf(stderr, "%.1fs\n", exec_time);
  //print_result();
  free(par);
  
  return 0;
}
