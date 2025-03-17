
#ifndef __UTILS_H
#define __UTILS_H

#include <stdbool.h>

bool cast_to_long(const char *str, long *out);
bool cast_to_long_long(const char *str, long long *out);
bool cast_to_double(const char *str, double *out);

#endif // __UTILS_H