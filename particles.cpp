
#include "particles.h"
#include "init_particles.h"
#include "math.h"

long find_particle_cell(const particle_t & par, double size, long ncside) {
  long xpart = par.x / size;
  long ypart = par.y / size;
  return ypart * ncside + xpart;
}

void copy_particle(particle_t & par1, const particle_t & par2) {
  par1.x = par2.x;
  par1.y = par2.y;
  par1.vx = par2.vx;
  par1.vy = par2.vy;
  par1.m = par2.m;
  par1.ind = par2.ind;
}

void remove_and_swap(std::vector<cell_t>& cells, long i, long long j) {
  long long last = cells[i].par.size() - 1;
  copy_particle(cells[i].par[j], cells[i].par[last]);
  cells[i].par.resize(last);
}

double calculate_distance(const particle_t & par1, const particle_t & par2) {
  double dx = par1.x - par2.x;
  double dy = par1.y - par2.y;
  return dx * dx + dy * dy;
}

void update_position_and_velocity(particle_t& par, const vec_t& acc) {
  par.x += par.vx + 0.5*(DELTAT*DELTAT)*acc.x;
  par.y += par.vy + 0.5*(DELTAT*DELTAT)*acc.y;
  par.vx += DELTAT*acc.x;
  par.vy += DELTAT*acc.y;
}

vec_t compute_gravitacional_force(const particle_t & par1, const particle_t & par2) {
  vec_t force;
  double dx = par1.x - par2.x;
  double dy = par1.y - par2.y;
  double distance = dx * dx + dy * dy;
  double hypotenuse = sqrt(distance);
  double cos = dx/hypotenuse;
  double sin = dy/hypotenuse;
  force.x = cos*(G * par1.m * par2.m) / distance;
  force.y = sin*(G * par1.m * par2.m) / distance;
  return force;
}
