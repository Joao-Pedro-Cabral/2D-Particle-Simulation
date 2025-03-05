
#include "particles.h"
#include "init_particles.h"
#include "math.h"

long find_particle_cell(const particle_t &par, double size, long ncside) {
  long xpart = par.x / size;
  long ypart = par.y / size;
  return INDEX(xpart, ypart, ncside);
}

void copy_particle(particle_t &par1, const particle_t &par2) {
  par1.x = par2.x;
  par1.y = par2.y;
  par1.vx = par2.vx;
  par1.vy = par2.vy;
  par1.ax = par2.ax;
  par1.ay = par2.ay;
  par1.m = par2.m;
  par1.ind = par2.ind;
  par1.collided = par2.collided;
}

double calc_squared_distance(const particle_t &par1, const particle_t &par2) {
  double dx = par1.x - par2.x;
  double dy = par1.y - par2.y;
  return dx * dx + dy * dy;
}

void update_position_and_velocity(double side, particle_t &par) {
  par.x += DELTAT * par.vx + 0.5 * (DELTAT * DELTAT) * par.ax;
  if (par.x > side)
    par.x -= side;
  if (par.x < 0)
    par.x += side;
  par.y += DELTAT * par.vy + 0.5 * (DELTAT * DELTAT) * par.ay;
  if (par.y > side)
    par.y -= side;
  if (par.y < 0)
    par.y += side;
  par.vx += DELTAT * par.ax;
  par.vy += DELTAT * par.ay;
  par.ax = 0.0;
  par.ay = 0.0;
}

vec_t calc_gravitational_force(const particle_t &par1, const particle_t &par2) {
  vec_t force;
  double dx = par2.x - par1.x;
  double dy = par2.y - par1.y;
  double denominator = dx * dx + dy * dy;
  denominator *= sqrt(denominator);
  double numerator = G * par1.m * par2.m;
  double F = numerator;
  force.x = dx * F;
  force.y = dy * F;
  return force;
}
