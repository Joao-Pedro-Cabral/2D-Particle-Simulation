
#ifndef __UTILS_H
#define __UTILS_H

#include <stdbool.h>
#include <cstdint>

bool cast_to_long(const char *str, uint32_t *out);
bool cast_to_long_long(const char *str, uint64_t *out);
bool cast_to_double(const char *str, double *out);

#endif // __UTILS_H