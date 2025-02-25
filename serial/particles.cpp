
#include "particles.h"
#include "init_particles.h"
#include "math.h"

uint32_t find_particle_cell(const particle_t &par, double size, uint32_t ncside) {
  uint32_t xpart = par.x / size;
  uint32_t ypart = par.y / size;
  return INDEX(xpart, ypart, ncside);
}

void copy_particle(particle_t &par1, const particle_t &par2) {
  par1.x = par2.x;
  par1.y = par2.y;
  par1.vx = par2.vx;
  par1.vy = par2.vy;
  par1.m = par2.m;
  par1.ind = par2.ind;
}

void remove_and_swap(std::vector<cell_t> &cells, uint32_t i, uint64_t j) {
  uint64_t last = cells[i].par.size() - 1;
  copy_particle(cells[i].par[j], cells[i].par[last]);
  cells[i].par.resize(last);
}

double calc_squared_distance(const particle_t &par1, const particle_t &par2) {
  double dx = par1.x - par2.x;
  double dy = par1.y - par2.y;
  return dx * dx + dy * dy;
}

void update_position_and_velocity(double side, particle_t &par,
                                  const vec_t &acc) {
  par.x += DELTAT * par.vx + 0.5 * (DELTAT * DELTAT) * acc.x;
  if (par.x > side)
    par.x -= side;
  if (par.x < 0)
    par.x += side;
  par.y += DELTAT * par.vy + 0.5 * (DELTAT * DELTAT) * acc.y;
  if (par.y > side) 
    par.y -= side;
  if (par.y < 0)
    par.y += side;
  par.vx += DELTAT * acc.x;
  par.vy += DELTAT * acc.y;
}

vec_t calc_gravitational_force(const particle_t &par1, const particle_t &par2) {
  vec_t force;
  double dx = par2.x - par1.x;
  double dy = par2.y - par1.y;
  double squared_distance = dx * dx + dy * dy;
  double hypotenuse = sqrt(squared_distance);
  double cos = dx / hypotenuse;
  double sin = dy / hypotenuse;
  force.x = cos * (G * par1.m * par2.m) / squared_distance;
  force.y = sin * (G * par1.m * par2.m) / squared_distance;
  return force;
}
