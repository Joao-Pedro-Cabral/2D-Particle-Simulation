
#ifndef __DEBUG_H
#define __DEBUG_H
#ifdef DEBUG_H
    #define DEBUG(str, ...) printf(str, ##__VA_ARGS__)
#else
    #define DEBUG(str, ...) ((void)0)
#endif // DEBUG_H
#endif // __DEBUG_H
