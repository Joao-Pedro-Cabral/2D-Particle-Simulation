
#include "init_particles.h"
#include "nodes.h"
#include "debug.h"
#include <math.h>

unsigned int seed;
void init_r4uni(int input_seed) { seed = input_seed + 987654321; }
double rnd_uniform01() {
  int seed_in = seed;
  seed ^= (seed << 13);
  seed ^= (seed >> 17);
  seed ^= (seed << 5);
  return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}
double rnd_normal01() {
  double u1, u2, z, result;
  do {
    u1 = rnd_uniform01();
    u2 = rnd_uniform01();
    z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
  } while (result < 0 || result >= 1);
  return result;
}

void init_particles(long seed, double side, long ncside, long long n_part,
                    particles_buffer_t *par, communication_buffers_t * buffers) {
  double (*rnd01)() = rnd_uniform01;
  long long i;

  if (seed < 0) {
    rnd01 = rnd_normal01;
    seed = -seed;
  }

  init_r4uni(seed);
  double size = side / ncside;

  for (i = 0; i < n_part; i++) {
    double x = rnd01();
    double y = rnd01();
    double vx = rnd01();
    double vy = rnd01();
    double m = rnd01();
    x = x * side;
    y = y * side;
    long xpart = x / size;
    long ypart = y / size;
    int col = (buffers->dims[1]*(xpart+1)-1) / ncside;
    int row = (buffers->dims[0]*(ypart+1)-1) / ncside;
    if(buffers->coords[0] == row && buffers->coords[1] == col) {
      vx = (vx - 0.5) * side / ncside / 5.0;
      vy = (vy - 0.5) * side / ncside / 5.0;
      m = m * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
      particles_buffer_add_p(par, x, y, vx, vy, m, i);
    }
  }
}
