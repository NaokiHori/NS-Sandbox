#include <stddef.h> // size_t
#include "logger.h" // LOGGER_FAILURE
#include "domain.h" // X_PERIODIC, Y_PERIODIC, nx, ny
#include "exchange_halo.h"

int exchange_halo_x(
    const domain_t * const domain,
    double ** const array
) {
  if (!X_PERIODIC) {
    LOGGER_FAILURE("x direction is not periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 0; j <= ny + 1; j++) {
    array[j][     0] = array[j][nx];
    array[j][nx + 1] = array[j][ 1];
  }
  return 0;
abort:
  return 1;
}

int exchange_halo_y(
    const domain_t * const domain,
    double ** const array
) {
  if (!Y_PERIODIC) {
    LOGGER_FAILURE("y direction is not periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t i = 0; i <= nx + 1; i++) {
    array[     0][i] = array[ny][i];
    array[ny + 1][i] = array[ 1][i];
  }
  return 0;
abort:
  return 1;
}

