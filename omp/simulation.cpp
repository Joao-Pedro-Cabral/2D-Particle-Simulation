#include "simulation.h"
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

void debug_center(long ncside, const std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside + 2; iy++) {
    for (long ix = 0; ix < ncside + 2; ix++) {
      long i = iy * (ncside + 2) + ix;
      DEBUG("Cell %ld x: %.6lf y: %.6lf m: %.6lf\n", i, cells[i].center.x,
            cells[i].center.y, cells[i].center.m);
    }
  }
}

void compute_centers_of_mass(double side, long ncside,
                             std::vector<cell_t> &cells) {
  #pragma omp parallel for collapse(2) //TODO ask professor if it is better to use a separate clause for parallel/for
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      double total_mass = 0.0;
      double weighted_x = 0.0;
      double weighted_y = 0.0;

      #pragma omp for reduction (+: total_mass, weighted_x, weighted_y)
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
              "mass to (-2side, -2side)\n",
              i);
        cells[i].center.x = -2 * side;
        cells[i].center.y = -2 * side;
      }
    }
  }

  #pragma omp parallel for collapse(2) nowait
  for (long i = 0; i < ncside + 2; i += ncside + 1) {
    for (long j = 0; j < ncside + 2; j++) {
      long ind = i * (ncside + 2) + j;
      long i2 = (i == 0) ? ncside : 1;
      if (j == 0 || j == ncside + 1) {
        long j2 = (j == 0) ? ncside : 1;
        long ind2 = i2 * (ncside + 2) + j2;
        cells[ind].center.x = (j == 0) ? cells[ind2].center.x - side
                                       : cells[ind2].center.x + side;
        cells[ind].center.y = (i == 0) ? cells[ind2].center.y - side
                                       : cells[ind2].center.y + side;
        cells[ind].center.m = cells[ind2].center.m;
      } else {
        long ind2 = i2 * (ncside + 2) + j;
        cells[ind].center.x = cells[ind2].center.x;
        cells[ind].center.y = (i == 0) ? cells[ind2].center.y - side
                                       : cells[ind2].center.y + side;
        cells[ind].center.m = cells[ind2].center.m;
      }
    }
  }

  #pragma omp parallel for collapse(2) 
  for (long j = 0; j < ncside + 2; j += ncside + 1) {
    for (long i = 1; i < ncside + 1; i++) {
      long ind = i * (ncside + 2) + j;
      long j2 = (j == 0) ? ncside : 1;
      long ind2 = i * (ncside + 2) + j2;
      cells[ind].center.x =
          (j == 0) ? cells[ind2].center.x - side : cells[ind2].center.x + side;
      cells[ind].center.y = cells[ind2].center.y;
      cells[ind].center.m = cells[ind2].center.m;
    }
  }
  debug_center(ncside, cells);
}

void debug_accelerations(const vec_t &force, long ind) {
  DEBUG("Force x: %.6lf, y: %.6lf, ind : %ld\n", force.x, force.y, ind);
}

void compute_accelerations(long ncside, std::vector<cell_t> &cells) {
  #pragma omp parallel for collapse(3)
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      for (long long j = 0; j < cells[INDEX(ix, iy, ncside)].par.size(); j++) {
        long i = INDEX(ix, iy, ncside);
        double resultant_x = 0.0;
        double resultant_y = 0.0;

        #pragma omp for reduction (+: resultant_x, resultant_y)
        for (long long k = 0; k < cells[i].par.size(); k++) {
          if (k == j)
            continue;
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[i].par[k]);
          resultant_x += force.x;
          resultant_y += force.y;
          if (cells[i].par[j].ind == 0)
            debug_accelerations(force, 0);
        }

        #pragma omp for reduction (+: resultant_x, resultant_y) // TOTODODO test wether performs better or not 
        for (long long k = 0; k < 9; k++) {
          if (k == 4)
            continue;
          long ind = i + ((k % 3) - 1) + (k / 3 - 1) * (ncside + 2);
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[ind].center);
          resultant_x += force.x;
          resultant_y += force.y;
          if (cells[i].par[j].ind == 0)
            debug_accelerations(force, ind);
        }
        cells[i].par[j].ax = resultant_x / cells[i].par[j].m;
        cells[i].par[j].ay = resultant_y / cells[i].par[j].m;
      }
    }
  }
}

void debug_compute_new_positions_and_velocities(long ncside,
                                                std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Particle %lld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: %.6f, "
              "cell: %ld\n",
              cells[i].par[j].ind, cells[i].par[j].m, cells[i].par[j].x,
              cells[i].par[j].y, cells[i].par[j].vx, cells[i].par[j].vy,
              iy * ncside + ix);
      }
    }
  }
}

void compute_new_positions_and_velocities(double side, double size, long ncside,
                                          std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        update_position_and_velocity(side, cells[i].par[j]);
      }
    }
  }

  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        long ind = find_particle_cell(cells[i].par[j], size, ncside);
        if (ind == i)
          continue;
        cells[ind].par.push_back(cells[i].par[j]);
        remove_and_swap(cells, i, j);
      }
    }
  }
  debug_compute_new_positions_and_velocities(ncside, cells);
}

long long detect_collisions(long ncside, std::vector<cell_t> &cells,
                            std::vector<bool> &collisions) {
  long long n_collisions = 0;
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Size: %ld, i: %ld, j: %lld\n", cells[i].par.size(),
              ix + iy * ncside, cells[i].par[j].ind);
        collisions[j] = false;
        for (long long k = 0; k < j; k++) {
          double squared_distance =
              calc_squared_distance(cells[i].par[j], cells[i].par[k]);
          DEBUG("Distance: %.6lf, ix: %ld, iy: %ld, j: %lld, k: %lld\n",
                squared_distance, ix, iy, cells[i].par[j].ind,
                cells[i].par[k].ind);
          if (squared_distance > EPSILON2)
            continue;
          if (collisions[k] == false) {
            n_collisions++;
            collisions[k] = true;
          }
          collisions[j] = true;
        }
      }
      for (long long j = 0; j < cells[i].par.size(); j++) {
        if (collisions[j] == false)
          continue;
        DEBUG("Removing Particle %lld from cells[%ld].par[%lld]\n",
              cells[i].par[j].ind, i, j);
        remove_and_swap(cells, i, j);
      }
    }
  }
  DEBUG("Number of collisions %lld\n", n_collisions);
  return n_collisions;
}

particle_t find_particle_zero(long ncside, std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        if (cells[i].par[j].ind == 0)
          return cells[i].par[j];
      }
    }
  }
  ERROR("Particle zero not found\n");
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep,
                             const std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells((ncside + 2) * (ncside + 2));
  std::vector<vec_t> accs(npart);
  std::vector<bool> collisions(npart);
  simulation_result res;
  res.number_of_collisions = 0;
  DEBUG("Size: %lf", size);
  fill_cells(size, ncside, npart, par, cells);
  for (long long i = 0; i < nstep; i++) {
    DEBUG("--------STEP: %lld --------------\n", i);
    compute_centers_of_mass(side, ncside, cells);
    compute_accelerations(ncside, cells, accs);
    compute_new_positions_and_velocities(side, size, ncside, cells, accs);
    res.number_of_collisions += detect_collisions(ncside, cells, collisions);
  }
  res.particle_zero = find_particle_zero(ncside, cells);
  return res;
}
