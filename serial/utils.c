
#include "utils.h"
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>


bool cast_to_long(const char *str, uint32_t *out) {
  if (str == NULL || *str == '\0')
    return false;

  char *endptr = NULL;
  errno = 0;
  uint32_t val = strtol(str, &endptr, 10);

  if (errno == ERANGE)
    return false;

  *out = val;
  return true;
}

bool cast_to_long_long(const char *str, uint64_t *out) {
  if (str == NULL || *str == '\0')
    return false;

  char *endptr = NULL;
  errno = 0;
  uint64_t val = strtoll(str, &endptr, 10);

  if (errno == ERANGE)
    return false;

  *out = val;
  return true;
}

bool cast_to_double(const char *str, double *out) {
  if (str == NULL || *str == '\0')
    return false;

  char *endptr = NULL;
  errno = 0;
  double val = strtod(str, &endptr);
  if (errno == ERANGE)
    return false;

  *out = val;
  return true;
}
