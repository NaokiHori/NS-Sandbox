#include "logger.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "impose_bc.h"
#include "exchange_halo.h"
#include "./internal.h"

int update_pressure (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
) {
  const array_t * const psi = flow_solver->psi;
  array_t * const p = flow_field->p;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      p[j][i] += psi[j][i];
    }
  }
  // exchange halo / impose bc
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(p)) {
      LOGGER_FAILURE("failed to exchange halo in x");
      goto abort;
    }
  } else {
    if (0 != impose_bc_p_x(p)) {
      LOGGER_FAILURE("failed to impose boundary condition in x");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(p)) {
      LOGGER_FAILURE("failed to exchange halo in y");
      goto abort;
    }
  } else {
    if (0 != impose_bc_p_y(p)) {
      LOGGER_FAILURE("failed to impose boundary condition in y");
      goto abort;
    }
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to update pressure field");
  return 1;
}

