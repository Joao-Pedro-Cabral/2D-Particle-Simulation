#include "particles.h"

void center_of_mass(double side, long ncside, long long n_part, particle_t *par, cell_t *cells){

    double size = side/ncside;


    for (long long i = 0; i < ncside*ncside; i++){
        cells[i].n_part = 0;
        cells[i].c_part = 0;
    }

    for (long long i = 0; i < n_part; i++){
        long xpart = par[i].x/size;
        long ypart = par[i].y/size;
        cells[ypart*ncside + xpart].n_part ++;
    }

    for (long long i = 0; i < ncside*ncside; i++){
        cells[i].par = malloc(cells[i].n_part*sizeof(particle_t));
    }

    for(long long i = 0; i < n_part; i++){
        long xpart = par[i].x/size;
        long ypart = par[i].y/size;
        long ind = ypart*ncside + xpart;
        cells[ind].par[cells[ind].c_part].x = par[i].x;
        cells[ind].par[cells[ind].c_part].y = par[i].y;
        cells[ind].par[cells[ind].c_part].vx = par[i].vx;
        cells[ind].par[cells[ind].c_part].vy = par[i].vy;
        cells[ind].par[cells[ind].c_part].m = par[i].m;
    }

}