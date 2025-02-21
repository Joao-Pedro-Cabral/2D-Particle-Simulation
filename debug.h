
#ifndef __DEBUG_H
#define __DEBUG_H
#ifdef DEBUG_H
#define DEBUG(str, ...) printf(str, ##__VA_ARGS__)
#define ERROR(str, ...)                                                        \
  do {                                                                         \
    printf(str, ##__VA_ARGS__);                                                \
    exit(1);                                                                   \
  } while (0)
#else
#define DEBUG(str, ...) ((void)0)
#define ERROR(str, ...) exit(1)
#endif // DEBUG_H
#endif // __DEBUG_H
