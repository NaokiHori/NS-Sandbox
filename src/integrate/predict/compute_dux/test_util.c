#if defined(TEST)

#include <assert.h>
#include <math.h>
#include "domain.h"
#include "./test_util.h"

static const double pi = 3.141592653589793;

double get_x (
    const size_t i
) {
  assert(0 < i);
  return 0.5 * (2 * i - 2) * DX;
}

double get_y (
    const size_t j
) {
  assert(0 < j);
  return 0.5 * (2 * j - 1) * DY;
}

double get_duxdx (
    const double x,
    const double y
) {
  return 2. * pi / LX * cos(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

double get_duxdy (
    const double x,
    const double y
) {
  return 2. * pi / LY * sin(2. * pi * x / LX) * cos(2. * pi * y / LY);
}

double get_d2uxdx2 (
    const double x,
    const double y
) {
  return - 4. * pi * pi / LX / LX * sin(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

double get_d2uxdy2 (
    const double x,
    const double y
) {
  return - 4. * pi * pi / LY / LY * sin(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

double get_dpdx (
    const double x,
    const double y
) {
  return 2. * pi / LX * cos(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

#else

extern char dummy;

#endif // TEST
