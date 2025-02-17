
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

bool cast_to_long(const char *str, long *out) {
  if (str == NULL || *str == '\0') return false;

  char *endptr;
  errno = 0;
  long val = strtol(str, &endptr, 10);

  if (errno == ERANGE) return false;

  *out = val;
  return true;

}

bool cast_to_long_long(const char *str, long long *out) {
  if (str == NULL || *str == '\0') return false;

  char *endptr;
  errno = 0;
  long long val = strtoll(str, &endptr, 10);

  if (errno == ERANGE) return false;

  *out = val;
  return true;

}

bool cast_to_double(const char *str, double *out) {
  if (str == NULL || *str == '\0') return false;

  char *endptr;
  errno = 0;
  double val = strtod(str, &endptr);
  if (errno == ERANGE) return false;

  *out = val;
  return true;
}

