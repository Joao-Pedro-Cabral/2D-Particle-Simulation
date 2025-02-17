#include <stdio.h>
#include "debug.h"
#include "utils.h"
#include "init_particles.h"

int main(int argc, char* argv[]) {
  if(argc != 5) {
    DEBUG("Incorrect number of arguments");
    return 1;
  }

  long seed;
  double side;
  long ncside;
  long long npart, nstep;

  if(!cast_to_long(argv[0], &seed)) {
    DEBUG("Unexpected format for seed");
    return 1;
  }

  if(!cast_to_double(argv[1], &side)) {
    DEBUG("Unexpected format for side");
    return 1;
  }

  if(!cast_to_long(argv[2], &ncside)) {
    DEBUG("Unexpected format for ncside");
    return 1;
  }

  if(!cast_to_long_long(argv[3], &npart)) {
    DEBUG("Unexpected format for npart");
    return 1;
  }

  if(!cast_to_long_long(argv[4], &nstep)) {
    DEBUG("Unexpected format for nstep");
    return 1;
  }


  return 0;
}