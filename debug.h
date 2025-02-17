
#ifndef __DEBUG_H
#define __DEBUG_H
#ifdef DEBUG_H
    #define DEBUG(str) printf("%s\n", str)
#else
    #define DEBUG(str) ((void)0)
#endif // DEBUG_H
#endif // __DEBUG_H
