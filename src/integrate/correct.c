#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "exchange_halo.h"
#include "./internal.h"

static int correct_ux (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const array_t * const psi = flow_solver->psi;
  array_t * const ux = flow_field->ux;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      ux[j][i] -= dt / DX * (
          - psi[j    ][i - 1]
          + psi[j    ][i    ]
      );
    }
  }
  // NOTE: since the scalar pressure does not modify velocities on the boundaries,
  //       only halo exchanges are done here (not imposing BCs again)
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(ux)) {
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(ux)) {
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int correct_uy (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const array_t * const psi = flow_solver->psi;
  array_t * const uy = flow_field->uy;
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      uy[j][i] -= dt / DY * (
          - psi[j - 1][i    ]
          + psi[j    ][i    ]
      );
    }
  }
  // NOTE: since the scalar pressure does not modify velocities on the boundaries,
  //       only halo exchanges are done here (not imposing BCs again)
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(uy)) {
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(uy)) {
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

int correct (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  if (0 != correct_ux(flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct uy");
    goto abort;
  }
  if (0 != correct_uy(flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct uy");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to correct velocity field");
  return 1;
}

