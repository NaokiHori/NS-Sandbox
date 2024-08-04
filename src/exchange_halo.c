#include <stddef.h> // size_t
#include "logger.h" // LOGGER_FAILURE
#include "domain.h" // X_PERIODIC, Y_PERIODIC, NX, NY
#include "exchange_halo.h"

int exchange_halo_x (
    array_t * const array
) {
  if (!X_PERIODIC) {
    LOGGER_FAILURE("x direction is not periodic");
    goto abort;
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    array[j][     0] = array[j][NX];
    array[j][NX + 1] = array[j][ 1];
  }
  return 0;
abort:
  return 1;
}

int exchange_halo_y (
    array_t * const array
) {
  if (!Y_PERIODIC) {
    LOGGER_FAILURE("y direction is not periodic");
    goto abort;
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    array[     0][i] = array[NY][i];
    array[NY + 1][i] = array[ 1][i];
  }
  return 0;
abort:
  return 1;
}

