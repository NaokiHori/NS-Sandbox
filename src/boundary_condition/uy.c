#include <stddef.h> // size_t
#include "logger.h"
#include "boundary_condition.h"

int impose_boundary_condition_uy_x(
    const domain_t * const domain,
    double ** const uy
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("x direction is periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double uy_xm = 0.;
  const double uy_xp = 0.;
  for (size_t j = 0; j <= ny + 1; j++) {
    uy[j][     0] = 2. * uy_xm - uy[j][ 1];
    uy[j][nx + 1] = 2. * uy_xp - uy[j][nx];
  }
  return 0;
abort:
  return 1;
}

int impose_boundary_condition_uy_y(
    const domain_t * const domain,
    double ** const uy
) {
  if (Y_PERIODIC) {
    LOGGER_FAILURE("y direction is periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t i = 0; i <= nx + 1; i++) {
    uy[     0][i] = uy[1][i];
    uy[ny + 1][i] = -1.;
  }
  return 0;
abort:
  return 1;
}

