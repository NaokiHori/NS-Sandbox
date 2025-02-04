#if defined(TEST)

#include <assert.h>
#include <math.h>
#include "./test_util.h"

static const double pi = 3.141592653589793;

double get_x(
    const domain_t * const domain,
    const size_t i
) {
  const double dx = domain->dx;
  assert(0 < i);
  return 0.5 * (2 * i - 1) * dx;
}

double get_y(
    const domain_t * const domain,
    const size_t j
) {
  const double dy = domain->dy;
  assert(0 < j);
  return 0.5 * (2 * j - 2) * dy;
}

double get_duydx(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return
    2. * pi / lx
    * cos(2. * pi * x / lx)
    * sin(2. * pi * y / ly);
}

double get_duydy(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return
    2. * pi / ly
    * sin(2. * pi * x / lx)
    * cos(2. * pi * y / ly);
}

double get_d2uydx2(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return
    - 4. * pi * pi / lx / lx
    * sin(2. * pi * x / lx)
    * sin(2. * pi * y / ly);
}

double get_d2uydy2(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return
    - 4. * pi * pi / ly / ly
    * sin(2. * pi * x / lx)
    * sin(2. * pi * y / ly);
}

double get_dpdy(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return
    2. * pi / ly
    * sin(2. * pi * x / lx)
    * cos(2. * pi * y / ly);
}

#else

extern char dummy;

#endif // TEST
