#include "particles.h"
#include <vector>

void center_of_mass(double side, long ncside, long long n_part, const std::vector<particle_t>& par,
                    std::vector<cell_t>& cells) {

  double size = side / ncside;

  for (long long i = 0; i < n_part; i++) {
    long xpart = par[i].x / size;
    long ypart = par[i].y / size;
    long ind = ypart * ncside + xpart;
    cells[ind].par.push_back(par[i]);
  }
}