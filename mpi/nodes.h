
#ifndef __NODES_H
#define __NODES_H

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p)) * n
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_NUM_OF_NEIGHBORS(id, p, n) (2)
#define BLOCK_FRONTIER(id, p, n) (n)
#define NUM_COLUMNS(id, p, n) (n)
#define NUM_ROWS(id, p, n) BLOCK_SIZE(id, p, n) / NUM_COLUMNS(id, p, n)
#define BLOCK_NEIGHBORHOOD(id, p, n) (2 * (n + 2) + 2 * NUM_ROWS(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n * n))
#define TAG_CENTER_UP 0
#define TAG_CENTER_DOWN 1
#define TAG_PARTICLE_UP 2
#define TAG_PARTICLE_DOWN 3

#endif // __NODES_H
