
#include "particles.h"
#include "init_particles.h"
#include "math.h"

long find_particle_cell(const particle_t &par, double size, long ncside) {
  long xpart = par.x / size;
  long ypart = par.y / size;
  return INDEX(xpart, ypart, ncside);
}

long find_particle_cell(const cell_t &cell, long long i, double size,
                        long ncside) {
  long xpart = cell.x[i] / size;
  long ypart = cell.y[i] / size;
  return INDEX(xpart, ypart, ncside);
}

void copy_particle(particle_t &par1, const cell_t &cell, long long i) {
  par1.x = cell.x[i];
  par1.y = cell.y[i];
  par1.vx = cell.vx[i];
  par1.vy = cell.vy[i];
  par1.m = cell.m[i];
}

void copy_particle(cell_t &cell, long long i, long long j) {
  cell.x[i] = cell.x[j];
  cell.y[i] = cell.y[j];
  cell.vx[i] = cell.vx[j];
  cell.vy[i] = cell.vy[j];
  cell.ax[i] = cell.ax[j];
  cell.ay[i] = cell.ay[j];
  cell.m[i] = cell.m[j];
  cell.ind[i] = cell.ind[j];
  cell.collided[i] = cell.collided[j];
}

double calc_squared_distance(const cell_t &cell, long long i, long long j) {
  double dx = cell.x[i] - cell.x[j];
  double dy = cell.y[i] - cell.y[j];
  return dx * dx + dy * dy;
}

void update_position_and_velocity(double side, cell_t &cell, long long i) {
  cell.x[i] += DELTAT * cell.vx[i] + 0.5 * (DELTAT * DELTAT) * cell.ax[i];
  if (cell.x[i] > side)
    cell.x[i] -= side;
  if (cell.x[i] < 0)
    cell.x[i] += side;
  cell.y[i] += DELTAT * cell.vy[i] + 0.5 * (DELTAT * DELTAT) * cell.ay[i];
  if (cell.y[i] > side)
    cell.y[i] -= side;
  if (cell.y[i] < 0)
    cell.y[i] += side;
  cell.vx[i] += DELTAT * cell.ax[i];
  cell.vy[i] += DELTAT * cell.ay[i];
  cell.ax[i] = 0.0;
  cell.ay[i] = 0.0;
}

vec_t calc_gravitational_force(const cell_t &cell, long long i, long long j) {
  vec_t force;
  double dx = cell.x[j] - cell.x[i];
  double dy = cell.y[j] - cell.y[i];
  double denominator = dx * dx + dy * dy;
  denominator *= sqrt(denominator);
  double numerator = G * cell.m[i] * cell.m[j];
  double F = numerator / denominator;
  force.x = dx * F;
  force.y = dy * F;
  return force;
}

vec_t calc_gravitational_force(const cell_t &cell, long long i,
                               const particle_t &center) {
  vec_t force;
  double dx = center.x - cell.x[i];
  double dy = center.y - cell.y[i];
  double denominator = dx * dx + dy * dy;
  denominator *= sqrt(denominator);
  double numerator = G * cell.m[i] * center.m;
  double F = numerator / denominator;
  force.x = dx * F;
  force.y = dy * F;
  return force;
}

void cell_t::push_back(const particle_t &par, long long i) {
  x.push_back(par.x);
  y.push_back(par.y);
  vx.push_back(par.vx);
  vy.push_back(par.vy);
  ax.push_back(0.0);
  ay.push_back(0.0);
  m.push_back(par.m);
  ind.push_back(i);
  collided.push_back(false);
}

void cell_t::push_back(const cell_t &cell, long long i) {
  x.push_back(cell.x[i]);
  y.push_back(cell.y[i]);
  vx.push_back(cell.vx[i]);
  vy.push_back(cell.vy[i]);
  ax.push_back(cell.ax[i]);
  ay.push_back(cell.ay[i]);
  m.push_back(cell.m[i]);
  ind.push_back(cell.ind[i]);
  collided.push_back(cell.collided[i]);
}

void cell_t::resize(long long size) {
  x.resize(size);
  y.resize(size);
  vx.resize(size);
  vy.resize(size);
  ax.resize(size);
  ay.resize(size);
  m.resize(size);
  ind.resize(size);
  collided.resize(size);
}
