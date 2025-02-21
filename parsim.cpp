#include "parsim.h"
#include "debug.h"
#include "init_particles.h"
#include "particles.h"
#include <vector>

void fill_cells(double size, long ncside, long long n_part,
                const std::vector<particle_t> &par,
                std::vector<cell_t> &cells) {
  for (long long i = 0; i < n_part; i++) {
    cells[find_particle_cell(par[i], size, ncside)].par.push_back(par[i]);
  }
}

void debug_center(const std::vector<cell_t> &cells) {
  for (long i = 0; i < cells.size(); i++) {
    DEBUG("Cell %ld x: %.3lf y: %.3lf m: %.3lf\n", i, cells[i].center.x,
          cells[i].center.y, cells[i].center.m);
  }
}

void compute_centers_of_mass(double size, long ncside,
                             std::vector<cell_t> &cells) {

  for (long i = 0; i < cells.size(); i++) {
    double total_mass = 0.0;
    double weighted_x = 0.0;
    double weighted_y = 0.0;

    for (long long j = 0; j < cells[i].par.size(); j++) {
      total_mass += cells[i].par[j].m;
      weighted_x += cells[i].par[j].m * cells[i].par[j].x;
      weighted_y += cells[i].par[j].m * cells[i].par[j].y;
    }

    cells[i].center.m = total_mass;

    if (total_mass > 0) {
      cells[i].center.x = weighted_x / total_mass;
      cells[i].center.y = weighted_y / total_mass;
    } else {
      DEBUG("Warning: Cell %ld has no particles. Setting default center of "
            "mass to (0.0, 0.0)\n",
            i);
      cells[i].center.x = (i % ncside + 0.5) * size;
      cells[i].center.y = (i / ncside + 0.5) * size;
    }
  }
  debug_center(cells);
}

void debug_accelerations(const vec_t& force) {
  DEBUG("Force x: %.3lf, y: %.3lf\n", force.x, force.y);
}

void compute_accelerations(long ncside, const std::vector<cell_t> &cells,
                           std::vector<vec_t> &accs) {
  long long l = 0;
  // TODO: WRAP!!!!!!!!!!!!!!!!!!!!!!!!!
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      double resultant_x = 0.0;
      double resultant_y = 0.0;
      for (long long k = 0; k < cells[i].par.size(); k++) {
        if (k == j)
          continue;
        vec_t force =
            calc_gravitacional_force(cells[i].par[j], cells[i].par[k]);
        resultant_x += force.x;
        resultant_y += force.y;
        if(cells[i].par[j].ind == 0) debug_accelerations(force);
      }
      for (long long k = 0; k < 9; k++) {
        if (k == 5)
          continue;
        long long ind =
            (i + ((k % 3) - 1) + (k / 3 - 1) * ncside) % cells.size();
        vec_t force =
            calc_gravitacional_force(cells[i].par[j], cells[ind].center);
        resultant_x += force.x;
        resultant_y += force.y;
        if(cells[i].par[j].ind == 0) debug_accelerations(force);
      }
      accs[l].x = resultant_x / cells[i].par[j].m;
      accs[l].y = resultant_y / cells[i].par[j].m;
      l++;
    }
  }
}

void compute_new_positions_and_velocities(double side, double size, long ncside,
                                          std::vector<cell_t> &cells,
                                          const std::vector<vec_t> &accs) {
  long long l = 0;
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      update_position_and_velocity(side, cells[i].par[j], accs[l]);
      l++;
    }
  }

  for (long i = 0; i < cells.size(); i++) {
    for (long long j = cells[i].par.size() - 1; j > 0; j--) {
      long ind = find_particle_cell(cells[i].par[j], size, ncside);
      if (ind == i)
        continue;
      cells[ind].par.push_back(cells[i].par[j]);
      remove_and_swap(cells, i, j);
    }
  }
}

long long detect_collisions(std::vector<cell_t> &cells,
                            std::vector<bool> &collisions) {
  long long n_collisions = 0;
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      collisions[j] = false;
      for (long long k = 0; k < j; k++) {
        double distance = calc_distance(cells[i].par[j], cells[i].par[k]);
        if (distance > DELTAT * DELTAT)
          continue;
        if (collisions[i] == false)
          n_collisions++;
        collisions[j] = true;
        collisions[i] = true;
      }
    }
    for (long long j = 0; j < cells[i].par.size(); j++) {
      if (collisions[j] == false)
        continue;
      remove_and_swap(cells, i, j);
    }
  }
  return n_collisions;
}

particle_t find_particle_zero(std::vector<cell_t> &cells) {
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      if (cells[i].par[j].ind == 0)
        return cells[i].par[j];
    }
  }
  ERROR("Particle zero not found\n");
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep,
                             const std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells(ncside * ncside);
  std::vector<vec_t> accs(npart);
  std::vector<bool> collisions(npart);
  simulation_result res;
  res.number_of_collisions = 0;
  fill_cells(size, ncside, npart, par, cells);
  for (long long i = 0; i < nstep; i++) {
    compute_centers_of_mass(size, ncside, cells);
    compute_accelerations(ncside, cells, accs);
    compute_new_positions_and_velocities(side, size, ncside, cells, accs);
    res.number_of_collisions += detect_collisions(cells, collisions);
  }
  res.particle_zero = find_particle_zero(cells);
  return res;
}
