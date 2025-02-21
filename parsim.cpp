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
      DEBUG("Warning: Cell  %d has no particles. Setting default center of "
            "mass to (0.0, 0.0)",
            i);
      cells[i].center.x = (i % ncside + 0.5) *
                          size; // Default if no particles exist in the cell
      cells[i].center.y = (i / ncside + 0.5) * size;
    }
  }
}

void compute_accelerations(long ncside, const std::vector<cell_t> &cells,
                           std::vector<vec_t> &acc_vec) {
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
            compute_gravitacional_force(cells[i].par[j], cells[i].par[k]);
        resultant_x += force.x;
        resultant_y += force.y;
      }
      for (long long k = 0; k < 9; k++) {
        if (k == 5)
          continue;
        long long ind = i + ((k % 3) - 1) + (k / 3 - 1) * ncside;
        vec_t force =
            compute_gravitacional_force(cells[i].par[j], cells[ind].center);
        resultant_x += force.x;
        resultant_y += force.y;
      }
      acc_vec[l].x = resultant_x / cells[i].par[j].m;
      acc_vec[l].y = resultant_y / cells[i].par[j].m;
      l++;
    }
  }
}

void compute_new_positions_and_velocities(double size, long ncside,
                                          std::vector<cell_t> &cells,
                                          const std::vector<vec_t> &acc_vec) {
  long long l = 0;
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      update_position_and_velocity(cells[i].par[j], acc_vec[l]);
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

long long detect_collisions(std::vector<cell_t> &cells) {
  long long n_collisions = 0;
  for (long i = 0; i < cells.size(); i++) {
    for (long long j = 0; j < cells[i].par.size(); j++) {
      for (long long k = 0; k < cells[i].par.size(); k++) {
        if (k == j)
          continue;
        double distance = calculate_distance(cells[i].par[j], cells[i].par[k]);
        if (distance > DELTAT * DELTAT)
          continue;
        remove_and_swap(cells, i, j);
        remove_and_swap(cells, i, k);
        n_collisions++;
      }
    }
  }
  return n_collisions;
}

void simulation(double side, long ncside, long long npart, long long nstep,
                const std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells(ncside);
  std::vector<vec_t> acc_vec(npart);
  fill_cells(size, ncside, npart, par, cells);
  for (long long i = 0; i < nstep; i++) {
    compute_centers_of_mass(size, ncside, cells);
    compute_accelerations(ncside, cells, acc_vec);
    compute_new_positions_and_velocities(size, ncside, cells, acc_vec);
    detect_collisions(cells);
  }
}
