#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "impose_bc.h"

int impose_bc_uy_x (
    array_t * const uy
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("x direction is periodic");
    goto abort;
  }
  const double uy_xm = 0.;
  const double uy_xp = 0.;
  for (size_t j = 0; j <= NY + 1; j++) {
    uy[j][     0] = 2. * uy_xm + uy[j][ 1];
    uy[j][NX + 1] = 2. * uy_xp - uy[j][NX];
  }
  return 0;
abort:
  return 1;
}

int impose_bc_uy_y (
    array_t * const uy
) {
  if (Y_PERIODIC) {
    LOGGER_FAILURE("y direction is periodic");
    goto abort;
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    uy[     0][i] = 0.;
    uy[     1][i] = 0.;
    uy[NY + 1][i] = 0.;
  }
  return 0;
abort:
  return 1;
}

