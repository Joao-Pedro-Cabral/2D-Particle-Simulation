
#ifndef __NODES_H
#define __NODES_H

#define BLOCK_INDEX(ind, k, n) (ind + ((k % 3) - 1) + (k / 3 - 1) * n)
#define NUM_OF_NEIGHBORS 8
#define TAG_CENTER 0
#define TAG_PARTICLE NUM_OF_NEIGHBORS
#define TAG_PARTICLE_ZERO 2*NUM_OF_NEIGHBORS

#endif // __NODES_H
