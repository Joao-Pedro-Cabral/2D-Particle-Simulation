#include "simulation.h"
#include "debug.h"
#include "init_particles.h"
#include "particles.h"
#include <vector>

void fill_cells(double size, long ncside, long long n_part,
                const std::vector<particle_t> &par, std::vector<cell_t> &cells) {
  for (long long i = 0; i < n_part; i++) {
    cell_t & cell = cells[find_particle_cell(par[i], size, ncside)];
    cell.push_back(par[i], i);
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

  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      double total_mass = 0.0;
      double weighted_x = 0.0;
      double weighted_y = 0.0;
      cell_t & cell = cells[i];
      long long cell_size = cell.x.size();

      for (long long j = 0; j < cell_size; j++) {
        total_mass += cell.m[j];
        weighted_x += cell.m[j] * cell.x[j];
        weighted_y += cell.m[j] * cell.y[j];
      }

      cell.center.m = total_mass;

      if (total_mass > 0) {
        cell.center.x = weighted_x / total_mass;
        cell.center.y = weighted_y / total_mass;
      } else {
        cell.center.x = -2 * side;
        cell.center.y = -2 * side;
      }
    }
  }
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
}

void compute_kinetics(double side, long ncside, std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      cell_t & cell = cells[i];
      long long cell_size = cell.x.size();
      for (long long j = 0; j < cell_size; j++) {
        for (long long k = j + 1; k < cell_size; k++) {
          vec_t force = calc_gravitational_force(cell, j, k);
          cell.ax[j] += force.x;
          cell.ay[j] += force.y;
          cell.ax[k] -= force.x;
          cell.ay[k] -= force.y;
        }
        for (long long k = 0; k < 9; k++) {
          if (k == 4)
            continue;
          long ind = i + ((k % 3) - 1) + (k / 3 - 1) * (ncside + 2);
          vec_t force = calc_gravitational_force(cell, j, cells[ind].center);
          cell.ax[j] += force.x;
          cell.ay[j] += force.y;
        }
        cell.ax[j] /= cell.m[j];
        cell.ay[j] /= cell.m[j];
        update_position_and_velocity(side, cell, j);
      }
    }
  }
}

void debug_particles(long ncside, std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      long long cell_size = cells[i].x.size();
      for (long long j = 0; j < cell_size; j++) {
        DEBUG("Particle %lld: m: %.6f, x: %.6f, y: %.6f, vx: %.6f, vy: %.6f, "
              "ax: %.6f, ay: %.6f, collided: %d, cell: %ld\n",
              cells[i].ind[j], cells[i].m[j], cells[i].x[j],
              cells[i].y[j], cells[i].vx[j], cells[i].vy[j],
              cells[i].ax[j], cells[i].ay[j], (int) cells[i].collided[j],
              iy * ncside + ix);
      }
    }
  }
}

void compute_new_particle_cell(double size, long ncside,
                               std::vector<cell_t> &cells) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      cell_t & cell = cells[i];
      long long cell_size = cell.x.size();
      for (long long j = 0; j < cell_size; j++) {
        long ind = find_particle_cell(cell, j, size, ncside);
        if (ind == i)
          continue;
        cells[ind].push_back(cell, j);
        cell.m[j] = 0;
      }
    }
  }
}

void detect_collisions(long ncside, std::vector<cell_t> &cells,
                       long long &n_collisions) {
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      long long cell_collisions = 0;
      cell_t & cell = cells[i];
      long long cell_size = cell.x.size();
      for (long long j = cell_size - 1; j >= 0; j--) {
        if (cell.m[j] == 0) {
          DEBUG("Removing mass null Particle %lld from cells[%ld].par[%lld]\n",
                cell.ind[j], i, j);
          cell_size--;
          copy_particle(cell, cell_size, j);
        }
      }
      for (long long j = cell_size - 1; j >= 0; j--) {
        for (long long k = j - 1; k >= 0; k--) {
          double squared_distance = calc_squared_distance(cell, j, k);
          if (squared_distance > EPSILON2)
            continue;
          DEBUG("Distance: %.6lf, ix: %ld, iy: %ld, j: %lld, k: %lld\n",
                squared_distance, ix, iy, cell.ind[j],
                cell.ind[k]);
          if (cell.collided[k] == false) {
            cell_collisions++;
            cell.collided[k] = true;
          }
          cell.collided[j] = true;
        }
        if (cell.collided[j] == true) {
          DEBUG("Removing Particle %lld from cells[%ld].par[%lld]\n",
                cells[i].ind[j], i, j);
          cell_size--;
          copy_particle(cell, cell_size, j);
        }
      }
      n_collisions += cell_collisions;
      cell.resize(cell_size);
    }
  }
  DEBUG("Number of collisions %lld\n", n_collisions);
  debug_particles(ncside, cells);
}

particle_t find_particle_zero(long ncside, std::vector<cell_t> &cells) {
  particle_t par0;
  for (long iy = 0; iy < ncside; iy++) {
    for (long ix = 0; ix < ncside; ix++) {
      long i = INDEX(ix, iy, ncside);
      long long cell_size = cells[i].x.size();
      for (long long j = 0; j < cell_size; j++) {
        if (cells[i].ind[j] == 0) {
          copy_particle(par0, cells[i], j);
          return par0;
        }
      }
    }
  }
  ERROR("Particle 0 not found\n");
}

simulation_result simulation(double side, long ncside, long long npart,
                             long long nstep, std::vector<particle_t> &par) {
  double size = side / ncside;
  std::vector<cell_t> cells((ncside + 2) * (ncside + 2));
  simulation_result res;
  res.number_of_collisions = 0;
  double exec_time1 = 0.0, exec_time2 = 0.0, exec_time3 = 0.0, exec_time4 = 0.0;
  fill_cells(size, ncside, npart, par, cells);
  debug_particles(ncside, cells);
  for (long long i = 0; i < nstep; i++) {
    printf("--------STEP: %lld --------------\n", i);
    exec_time1 -= omp_get_wtime();
    compute_centers_of_mass(side, ncside, cells);
    exec_time1 += omp_get_wtime();
    exec_time2 -= omp_get_wtime();
    compute_kinetics(side, ncside, cells);
    exec_time2 += omp_get_wtime();
    exec_time3 -= omp_get_wtime();
    compute_new_particle_cell(size, ncside, cells);
    exec_time3 += omp_get_wtime();
    exec_time4 -= omp_get_wtime();
    detect_collisions(ncside, cells, res.number_of_collisions);
    exec_time4 += omp_get_wtime();
  }
  fprintf(stderr, "Centers: %.1fs\n", exec_time1);
  fprintf(stderr, "Kinetics: %.1fs\n", exec_time2);
  fprintf(stderr, "Find: %.1fs\n", exec_time3);
  fprintf(stderr, "Collisions: %.1fs\n", exec_time4);
  res.particle_zero = find_particle_zero(ncside, cells);
  return res;
}
