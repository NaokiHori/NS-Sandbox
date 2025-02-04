#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "boundary_condition.h"

int impose_boundary_condition_ux_x(
    array_t * const ux
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("x direction is periodic");
    goto abort;
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    ux[j][     0] = 0.;
    ux[j][     1] = 0.;
    ux[j][NX + 1] = 0.;
  }
  return 0;
abort:
  return 1;
}

int impose_boundary_condition_ux_y(
    array_t * const ux
) {
  if (Y_PERIODIC) {
    LOGGER_FAILURE("y direction is periodic");
    goto abort;
  }
  const double ux_ym = 0.;
  const double ux_yp = 0.;
  for (size_t i = 0; i <= NX + 1; i++) {
    ux[     0][i] = 2. * ux_ym + ux[ 1][i];
    ux[NY + 1][i] = 2. * ux_yp - ux[NY][i];
  }
  return 0;
abort:
  return 1;
}

