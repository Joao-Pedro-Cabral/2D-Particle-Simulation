
#ifndef __NODES_H
#define __NODES_H

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))*n
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_NEIGHBORHOOD(id,p,n) (4 * (n + 1))
#define BLOCK_NUM_OF_NEIGHBORS(id,p,n) (2)
#define BLOCK_FRONTIER(id,p,n) (n)
#define NUM_ROWS(id,p,n) ((id)*(n)/(p))
#define NUM_COLUMNS(id,p,n) (n)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n*n))
#define TAG_CENTER_TOP 0
#define TAG_CENTER_BOTTOM 1

#endif // __NODES_H
