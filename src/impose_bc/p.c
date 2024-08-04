#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "impose_bc.h"

int impose_bc_p_x (
    array_t * const p
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("x direction is periodic");
    goto abort;
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    p[j][     0] = p[j][ 1];
    p[j][NX + 1] = p[j][NX];
  }
  return 0;
abort:
  return 1;
}

int impose_bc_p_y (
    array_t * const p
) {
  if (Y_PERIODIC) {
    LOGGER_FAILURE("y direction is periodic");
    goto abort;
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    p[     0][i] = p[ 1][i];
    p[NY + 1][i] = p[NY][i];
  }
  return 0;
abort:
  return 1;
}

