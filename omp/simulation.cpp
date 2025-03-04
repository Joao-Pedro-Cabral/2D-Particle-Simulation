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

void init_cells_lock(long ncside, std::vector<cell_t> &cells) {
   for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      omp_init_lock(&cells[i].lock);
    }
  }
}

void destroy_cells_lock(long ncside, std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      omp_destroy_lock(&cells[i].lock);
    }
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
  #pragma omp for collapse(2) schedule(dynamic, CHUNK_SIZE)
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
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
        cells[i].center.x = -2 * side;
        cells[i].center.y = -2 * side;
      }
    }
  }

  #pragma omp for collapse(2) nowait
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

  #pragma omp for collapse(2) 
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
  // #pragma omp single
  // debug_center(ncside, cells);
}

void compute_kinetics(double side, long ncside, std::vector<cell_t> &cells) {
  #pragma omp for collapse(2) schedule(dynamic, CHUNK_SIZE)
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        double resultant_x = 0.0;
        double resultant_y = 0.0;

        for (long long k = 0; k < cells[i].par.size(); k++) {
          if (k == j)
            continue;
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[i].par[k]);
          resultant_x += force.x;
          resultant_y += force.y;
        }

        for (long long k = 0; k < 9; k++) {
          if (k == 4)
            continue;
          long ind = i + ((k % 3) - 1) + (k / 3 - 1) * (ncside + 2);
          vec_t force =
              calc_gravitational_force(cells[i].par[j], cells[ind].center);
          resultant_x += force.x;
          resultant_y += force.y;
        }
        cells[i].par[j].ax = resultant_x / cells[i].par[j].m;
        cells[i].par[j].ay = resultant_y / cells[i].par[j].m;
      }

      for (long long j = 0; j < cells[i].par.size(); j++) {
        update_position_and_velocity(side, cells[i].par[j]);
      }
    }
  }
}

void debug_particles(long ncside, std::vector<cell_t> &cells) {
  DEBUG("PARTICLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      for (long long j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Particle %lld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: %.6f, "
              "ax: %.6f, ay: %.6f, collided: %d, cell: %ld\n",
              cells[i].par[j].ind, cells[i].par[j].m, cells[i].par[j].x,
              cells[i].par[j].y, cells[i].par[j].vx, cells[i].par[j].vy,
              cells[i].par[j].ax, cells[i].par[j].ay, cells[i].par[j].collided,
              iy * ncside + ix);
      }
    }
  }
}

void compute_new_particle_cell(double size, long ncside,
                                          std::vector<cell_t> &cells) {
  #pragma omp for collapse(2) schedule(dynamic, CHUNK_SIZE)                                          
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      omp_set_lock(&cells[i].lock);
      long long cell_size = cells[i].par.size();
      omp_unset_lock(&cells[i].lock);
      for (long long j = 0; j < cell_size; j++) {
        long ind = find_particle_cell(cells[i].par[j], size, ncside);
        if (ind == i)
          continue;
        omp_set_lock(&cells[ind].lock);
        cells[ind].par.push_back(cells[i].par[j]);
        omp_unset_lock(&cells[ind].lock);
        cells[i].par[j].m = 0;
      }
    }
  }
  // #pragma omp single
  // debug_particles(ncside, cells);
}

void detect_collisions(long ncside, std::vector<cell_t> &cells, long long &n_collisions) {
  #pragma omp for collapse(2) reduction (+:n_collisions) schedule(dynamic, CHUNK_SIZE)
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      long long cell_collisions = 0;
      for (long long j = 0; j < cells[i].par.size(); j++) {
        if (cells[i].par[j].m == 0) {
          DEBUG("Removing mass null Particle %lld from cells[%ld].par[%lld]\n",
            cells[i].par[j].ind, i, j);
          remove_and_swap(cells, i, j);
          j--;
        }
      }
      for (long long j = 0; j < cells[i].par.size(); j++) {
        DEBUG("Size: %ld, i: %ld, j: %lld\n", cells[i].par.size(),
              ix + iy * ncside, cells[i].par[j].ind);
        cells[i].par[j].collided = false;
        for (long long k = 0; k < j; k++) {
          double squared_distance =
              calc_squared_distance(cells[i].par[j], cells[i].par[k]);
          if (squared_distance > EPSILON2)
            continue;
          DEBUG("Distance: %.6lf, ix: %ld, iy: %ld, j: %lld, k: %lld\n",
                squared_distance, ix, iy, cells[i].par[j].ind,
                cells[i].par[k].ind);
          if (cells[i].par[k].collided == false) {
            cell_collisions++;
            cells[i].par[k].collided = true;
          }
          cells[i].par[j].collided = true;
        }
      }
      n_collisions += cell_collisions;
  
      for (long long j = 0; j < cells[i].par.size(); j++) {
        if (cells[i].par[j].collided == false)
          continue;
        DEBUG("Removing Particle %lld from cells[%ld].par[%lld]\n",
              cells[i].par[j].ind, i, j);
        remove_and_swap(cells, i, j);
        j--;
      }
    }
  }
  // #pragma omp single
  // debug_particles(ncside, cells);
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
  ERROR("Particle 0 not found\n");
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep,
                             const std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells((ncside + 2) * (ncside + 2));
  simulation_result res;
  res.number_of_collisions = 0;
  fill_cells(size, ncside, npart, par, cells);
  init_cells_lock(ncside, cells);
  #pragma omp parallel 
  for (long long i = 0; i < nstep; i++) {
    // #pragma omp single
    // DEBUG("--------STEP: %lld --------------\n", i);
    compute_centers_of_mass(side, ncside, cells);
    compute_kinetics(side, ncside, cells);
    compute_new_particle_cell(size, ncside, cells);
    detect_collisions(ncside, cells, res.number_of_collisions);
  }
  destroy_cells_lock(ncside, cells);
  res.particle_zero = find_particle_zero(ncside, cells);
  return res;
}
