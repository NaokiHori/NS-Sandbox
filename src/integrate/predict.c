#include "logger.h"
#include "param.h"
#include "boundary_condition.h"
#include "exchange_halo.h"
#include "./predict.h"
#include "./predict/compute_dux.h"
#include "./predict/compute_duy.h"

static int update_ux(
    const domain_t * const domain,
    double ** const dux,
    double ** const weight,
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = ux_imin; i <= nx; i++) {
      ux[j][i] += dux[j][i];
    }
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = ux_imin; i <= nx; i++) {
      const double w_xm = weight[j    ][i - 1];
      const double w_xp = weight[j    ][i    ];
      ux[j][i] *= (
          + 0.5 * w_xm
          + 0.5 * w_xp
      );
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(domain, ux)) {
      LOGGER_FAILURE("failed to exchange halo in x (ux)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_ux_x(domain, ux)) {
      LOGGER_FAILURE("failed to impose boundary condition in x (ux)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(domain, ux)) {
      LOGGER_FAILURE("failed to exchange halo in y (ux)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_ux_y(domain, ux)) {
      LOGGER_FAILURE("failed to impose boundary condition in y (ux)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int update_uy(
    const domain_t * const domain,
    double ** const duy,
    double ** const weight,
    double ** const uy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      uy[j][i] += duy[j][i];
    }
  }
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double w_ym = weight[j - 1][i    ];
      const double w_yp = weight[j    ][i    ];
      uy[j][i] *= (
          + 0.5 * w_ym
          + 0.5 * w_yp
      );
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(domain, uy)) {
      LOGGER_FAILURE("failed to exchange halo in x (uy)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_uy_x(domain, uy)) {
      LOGGER_FAILURE("failed to impose boundary condition in x (uy)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(domain, uy)) {
      LOGGER_FAILURE("failed to exchange halo in y (uy)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_uy_y(domain, uy)) {
      LOGGER_FAILURE("failed to impose boundary condition in y (uy)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

int predict(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  double ** const dux = flow_solver->dux;
  double ** const duy = flow_solver->duy;
  if (0 != compute_dux(domain, flow_field, dt, dux)) {
    LOGGER_FAILURE("failed to find dux");
    goto abort;
  }
  if (0 != compute_duy(domain, flow_field, dt, duy)) {
    LOGGER_FAILURE("failed to find duy");
    goto abort;
  }
  if (0 != update_ux(domain, dux, flow_field->weight, flow_field->ux)) {
    LOGGER_FAILURE("failed to update ux");
    goto abort;
  }
  if (0 != update_uy(domain, duy, flow_field->weight, flow_field->uy)) {
    LOGGER_FAILURE("failed to update uy");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to predict flow field");
  return 1;
}

