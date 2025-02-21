#include "particles.h"
#include "init_particles.h"
#include "debug.h"
#include <vector>
#include <cmath>


void particles_per_cell(double side, long ncside, long long n_part, const std::vector<particle_t>& par,
               std::vector<cell_t>& cells) {
  double size = side / ncside;

  for (long long i = 0; i < n_part; i++) {
    long xpart = par[i].x / size;
    long ypart = par[i].y / size;
    long ind = ypart * ncside + xpart;
    cells[ind].par.push_back(par[i]);
  }
}

void center_of_mass(double side, long ncside, long long n_part, std::vector<cell_t>& cells) {

  for (long long i = 0; i <cells.size(); i++) {
    double total_mass = 0.0;
    double weighted_x = 0.0;
    double weighted_y = 0.0;

    for (long long j = 0; j < cells[i].par.size(); j++){
      total_mass += cells[i].par[j].m;
      weighted_x += cells[i].par[j].m * cells[i].par[j].x;
      weighted_y += cells[i].par[j].m * cells[i].par[j].y;
    }

    cells[i].m = total_mass;

    if (total_mass > 0) {
      cells[i].x = weighted_x / total_mass;
      cells[i].y = weighted_y / total_mass;
    } else {
      DEBUG("Warning: Cell  %d has no particles. Setting default center of mass to (0.0, 0.0)", i);
      cells[i].x = (i%ncside + 0.5)*side;  // Default if no particles exist in the cell
      cells[i].y = (i/ncside + 0.5)*side;
    }
  }
}

void gravitational_force(double side, long ncside, long long n_part, std::vector<cell_t>& cells, std::vector<acc_t>& acc_vec) {
  long long l = 0;
  for (long long i = 0; i < cells.size(); i++) {
    for (long long  j = 0; j < cells[i].par.size(); j++) {
      double force_x = 0.0;
      double force_y = 0.0;
      for(long long k = 0; k < cells[i].par.size(); k ++) {
        if(k ==  j) continue;
        double dx = cells[i].par[j].x - cells[i].par[k].x;
        double dy = cells[i].par[j].y - cells[i].par[k].y;
        double distance = dx * dx + dy * dy;
        double hypotenuse = sqrt(distance);
        double cos = dx/hypotenuse;
        double sin = dy/hypotenuse;
        force_x += cos*(G * cells[i].par[k].m * cells[i].par[j].m) / distance;
        force_y += sin*(G * cells[i].par[k].m * cells[i].par[j].m) / distance;
      }
      for(long long k = 0; k < 9; k ++) {
        if(k ==  5) continue;
        long long ind = i + ((k%3)-1) + (k/3-1)*ncside;
        double dx = cells[i].par[j].x - cells[ind].x;
        double dy = cells[i].par[j].y - cells[ind].y;
        double distance = dx * dx + dy * dy;
        double hypotenuse = sqrt(distance);
        double cos = dx/hypotenuse;
        double sin = dy/hypotenuse;
        force_x += cos*(G * cells[ind].m * cells[i].par[j].m) / distance;
        force_y += sin*(G * cells[ind].m * cells[i].par[j].m) / distance;
      }
      acc_vec[l].x = force_x/cells[i].par[j].m;
      acc_vec[l].y = force_y/cells[i].par[j].m;
      l ++;
    }
  }

  l = 0;
  for (long long i = 0; i < cells.size(); i++) {
    for (long long  j = 0; j < cells[i].par.size(); j++) {
      cells[i].par[j].x += cells[i].par[j].vx + 0.5*(DELTAT*DELTAT)*acc_vec[l].x;
      cells[i].par[j].y += cells[i].par[j].vy + 0.5*(DELTAT*DELTAT)*acc_vec[l].y;
      cells[i].par[j].vx += DELTAT*acc_vec[l].x;
      cells[i].par[j].vy += DELTAT*acc_vec[l].y;
      l ++;
    }
  }
  
}