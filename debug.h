
#ifndef __DEBUG_H
#define __DEBUG_H
#ifdef DEBUG_H
#define DEBUG(str, ...) printf(str, ##__VA_ARGS__)
#define ERROR(fmt, ...)                                                        \
  do {                                                                         \
    printf("ERROR: " fmt "\n", ##__VA_ARGS__);                                 \
    exit(1);                                                                   \
  } while (0)
#else
#define DEBUG(str, ...) ((void)0)
#define ERROR(fmt, ...) exit(1)
#endif // DEBUG_H
#endif // __DEBUG_H
