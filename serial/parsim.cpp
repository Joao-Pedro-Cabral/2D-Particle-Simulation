#include "parsim.h"
#include "debug.h"
#include "init_particles.h"
#include "particles.h"
#include <vector>

void fill_cells(double size, uint32_t ncside, uint64_t n_part,
                const std::vector<particle_t> &par,
                std::vector<cell_t> &cells) {
  for (uint64_t i = 0; i < n_part; i++) {
    cells[find_particle_cell(par[i], size, ncside)].par.push_back(par[i]);
  }
}

void debug_center(uint32_t ncside, const std::vector<cell_t> &cells) {
  for (uint32_t iy = 0; iy < ncside + 2; iy++) {
    for (uint32_t ix = 0; ix < ncside + 2; ix++) {
      uint32_t i = iy * (ncside + 2) + ix;
      DEBUG("Cell %d x: %.6lf y: %.6lf m: %.6lf\n", i, cells[i].center.x,
            cells[i].center.y, cells[i].center.m);
    }
  }
}

void compute_centers_of_mass(double side, uint32_t ncside,
                             std::vector<cell_t> &cells) {

  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      double total_mass = 0.0;
      double weighted_x = 0.0;
      double weighted_y = 0.0;

      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        total_mass += cells[i].par[j].m;
        weighted_x += cells[i].par[j].m * cells[i].par[j].x;
        weighted_y += cells[i].par[j].m * cells[i].par[j].y;
      }

      cells[i].center.m = total_mass;

      if (total_mass > 0) {
        cells[i].center.x = weighted_x / total_mass;
        cells[i].center.y = weighted_y / total_mass;
      } else {
        DEBUG("Warning: Cell %d has no particles. Setting default center of "
              "mass to (-2side, -2side)\n",
              i);
        cells[i].center.x = -2 * side;
        cells[i].center.y = -2 * side;
      }
    }
  }
  for (uint32_t i = 0; i < ncside + 2; i += ncside + 1) {
    for (uint32_t j = 0; j < ncside + 2; j++) {
      uint32_t ind = i * (ncside + 2) + j;
      uint32_t i2 = (i == 0) ? ncside : 1;
      if (j == 0 || j == ncside + 1) {
        uint32_t j2 = (j == 0) ? ncside : 1;
        uint32_t ind2 = i2 * (ncside + 2) + j2;
        cells[ind].center.x = (j == 0) ? cells[ind2].center.x - side
                                       : cells[ind2].center.x + side;
        cells[ind].center.y = (i == 0) ? cells[ind2].center.y - side
                                       : cells[ind2].center.y + side;
        cells[ind].center.m = cells[ind2].center.m;
      } else {
        uint32_t ind2 = i2 * (ncside + 2) + j;
        cells[ind].center.x = cells[ind2].center.x;
        cells[ind].center.y = (i == 0) ? cells[ind2].center.y - side
                                       : cells[ind2].center.y + side;
        cells[ind].center.m = cells[ind2].center.m;
      }
    }
  }
  for (uint32_t j = 0; j < ncside + 2; j += ncside + 1) {
    for (uint32_t i = 1; i < ncside + 1; i++) {
      uint32_t ind = i * (ncside + 2) + j;
      uint32_t j2 = (j == 0) ? ncside : 1;
      uint32_t ind2 = i * (ncside + 2) + j2;
      cells[ind].center.x =
          (j == 0) ? cells[ind2].center.x - side : cells[ind2].center.x + side;
      cells[ind].center.y = cells[ind2].center.y;
      cells[ind].center.m = cells[ind2].center.m;
    }
  }
  debug_center(ncside, cells);
}

void debug_accelerations(const vec_t &force, uint32_t ind) {
  DEBUG("Force x: %.6lf, y: %.6lf, ind : %d\n", force.x, force.y, ind);
}

void compute_accelerations(uint32_t ncside, const std::vector<cell_t> &cells,
                           std::vector<vec_t> &accs) {
  uint64_t l = 0;
  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        double resultant_x = 0.0;
        double resultant_y = 0.0;
        for (uint64_t k = 0; k < cells[i].par.size(); k++) {
          if (k == j)
            continue;
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[i].par[k]);
          resultant_x += force.x;
          resultant_y += force.y;
          if (cells[i].par[j].ind == 0)
            debug_accelerations(force, 0);
        }
        for (uint64_t k = 0; k < 9; k++) {
          if (k == 4)
            continue;
          uint32_t ind = i + ((k % 3) - 1) + (k / 3 - 1) * (ncside + 2);
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[ind].center);
          resultant_x += force.x;
          resultant_y += force.y;
          if (cells[i].par[j].ind == 0)
            debug_accelerations(force, ind);
        }
        accs[l].x = resultant_x / cells[i].par[j].m;
        accs[l].y = resultant_y / cells[i].par[j].m;
        l++;
      }
    }
  }
}

void debug_compute_new_positions_and_velocities(uint32_t ncside, std::vector<cell_t> &cells) {
  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Particle %ld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: %.6f, cell: %d\n", cells[i].par[j].ind, cells[i].par[j].m,
                                        cells[i].par[j].x, cells[i].par[j].y, cells[i].par[j].vx, cells[i].par[j].vy, iy*ncside + ix);
      }
    }
  }
}

void compute_new_positions_and_velocities(double side, double size, uint32_t ncside,
                                          std::vector<cell_t> &cells,
                                          const std::vector<vec_t> &accs) {
  uint64_t l = 0;
  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        update_position_and_velocity(side, cells[i].par[j], accs[l]);
        l++;
      }
    }
  }

  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        uint32_t ind = find_particle_cell(cells[i].par[j], size, ncside);
        if (ind == i)
          continue;
        cells[ind].par.push_back(cells[i].par[j]);
        remove_and_swap(cells, i, j);
      }
    }
  }
  debug_compute_new_positions_and_velocities(ncside, cells);
}

uint64_t detect_collisions(uint32_t ncside, std::vector<cell_t> &cells,
                            std::vector<bool> &collisions) {
  uint64_t n_collisions = 0;
  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Size: %ld, i: %d, j: %ld\n", cells[i].par.size(), ix + iy*ncside, cells[i].par[j].ind);
        collisions[j] = false;
        for (uint64_t k = 0; k < j; k++) {
          double distance = calc_squared_distance(cells[i].par[j], cells[i].par[k]);
          DEBUG("Distance: %.6lf, ix: %d, iy: %d, j: %ld, k: %ld\n", distance, ix, iy, cells[i].par[j].ind, cells[i].par[k].ind);
          if (distance > EPSILON2)
            continue;
          if (collisions[k] == false) {
            n_collisions++;
            collisions[k] = true;
          }
          collisions[j] = true;
        }
      }
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        if (collisions[j] == false)
          continue;
        DEBUG("Removing Particle %ld from cells[%d].par[%ld]\n", cells[i].par[j].ind, i, j);
        remove_and_swap(cells, i, j);
      }
    }
  }
  DEBUG("Number of collisions %ld\n", n_collisions);
  return n_collisions;
}

particle_t find_particle_zero(uint32_t ncside, std::vector<cell_t> &cells) {
  for (uint32_t iy = 0; iy < ncside; iy++) {
    for (uint32_t ix = 0; ix < ncside; ix++) {
      uint32_t i = INDEX(ix, iy, ncside);
      for (uint64_t j = 0; j < cells[i].par.size(); j++) {
        if (cells[i].par[j].ind == 0)
          return cells[i].par[j];
      }
    }
  }
  ERROR("Particle zero not found\n");
}

simulation_result simulation(double side, uint32_t ncside, uint64_t npart,
                             uint64_t nstep,
                             const std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells((ncside + 2) * (ncside + 2));
  std::vector<vec_t> accs(npart);
  std::vector<bool> collisions(npart);
  simulation_result res;
  res.number_of_collisions = 0;
  DEBUG("Size: %lf", size);
  fill_cells(size, ncside, npart, par, cells);
  for (uint64_t i = 0; i < nstep; i++) {
    DEBUG("--------STEP: %ld --------------\n", i);
    compute_centers_of_mass(side, ncside, cells);
    compute_accelerations(ncside, cells, accs);
    compute_new_positions_and_velocities(side, size, ncside, cells, accs);
    res.number_of_collisions += detect_collisions(ncside, cells, collisions);
  }
  res.particle_zero = find_particle_zero(ncside, cells);
  return res;
}
