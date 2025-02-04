#include "logger.h"
#include "exchange_halo.h"
#include "./update_pressure.h"

int update_pressure(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  double ** const psi = flow_solver->psi;
  double ** const p = flow_field->p;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      p[j][i] += psi[j][i];
    }
  }
  // exchange halo
  // NOTE: since DCT assumes dpdx = 0,
  //       boundary conditions are not directly imposed
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(domain, p)) {
      LOGGER_FAILURE("failed to exchange halo in x");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(domain, p)) {
      LOGGER_FAILURE("failed to exchange halo in y");
      goto abort;
    }
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to update pressure field");
  return 1;
}

