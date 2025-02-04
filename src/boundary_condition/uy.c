#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "boundary_condition.h"

int impose_boundary_condition_uy_x(
    array_t * const uy
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("x direction is periodic");
    goto abort;
  }
  const double uy_xm = 0.;
  const double uy_xp = 0.;
  for (size_t j = 0; j <= NY + 1; j++) {
    uy[j][     0] = 2. * uy_xm - uy[j][ 1];
    uy[j][NX + 1] = 2. * uy_xp - uy[j][NX];
  }
  return 0;
abort:
  return 1;
}

int impose_boundary_condition_uy_y(
    array_t * const uy
) {
  if (Y_PERIODIC) {
    LOGGER_FAILURE("y direction is periodic");
    goto abort;
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    uy[     0][i] = uy[1][i];
    uy[NY + 1][i] = -1.;
  }
  return 0;
abort:
  return 1;
}

