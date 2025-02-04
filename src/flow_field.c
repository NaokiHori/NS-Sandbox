#include <math.h>
#include "logger.h"
#include "array.h"
#include "domain.h"
#include "flow_field.h"
#include "boundary_condition.h"
#include "exchange_halo.h"

static int init_ux(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      ux[j][i] = 0.;
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

static int init_uy(
    const domain_t * const domain,
    double ** const uy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      uy[j][i] = -1.;
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

static int init_p(
    const domain_t * const domain,
    double ** const p
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      p[j][i] = 0.;
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(domain, p)) {
      LOGGER_FAILURE("failed to exchange halo in x (p)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(domain, p)) {
      LOGGER_FAILURE("failed to exchange halo in y (p)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_weight(
    const domain_t * const domain,
    double ** const weight
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double center[2] = {
    0.501 * lx,
    5. * ly / 6.,
  };
  const double radius = lx / 16.;
  for (size_t j = 0; j <= ny + 1; j++) {
    const double y = 0.5 * (2 * j - 1) * dy;
    for (size_t i = 0; i <= nx + 1; i++) {
      const double x = 0.5 * (2 * i - 1) * dx;
      const double d = sqrt(
          + pow(x - center[0], 2.)
          + pow(y - center[1], 2.)
      );
      weight[j][i] = 0.5 * (1. + tanh(ny * (d - radius)));
    }
  }
  return 0;
}

int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  array_init(nx + 2, ny + 2, &flow_field->ux);
  array_init(nx + 2, ny + 2, &flow_field->uy);
  array_init(nx + 2, ny + 2, &flow_field->p);
  array_init(nx + 2, ny + 2, &flow_field->weight);
  if (0 != init_ux(domain, flow_field->ux)) {
    LOGGER_FAILURE("failed to initialize ux");
    goto abort;
  }
  if (0 != init_uy(domain, flow_field->uy)) {
    LOGGER_FAILURE("failed to initialize uy");
    goto abort;
  }
  if (0 != init_p(domain, flow_field->p)) {
    LOGGER_FAILURE("failed to initialize p");
    goto abort;
  }
  if (0 != init_weight(domain, flow_field->weight)) {
    LOGGER_FAILURE("failed to initialize weight");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int flow_field_finalize(
    flow_field_t * const flow_field
) {
  array_finalize(&flow_field->ux);
  array_finalize(&flow_field->uy);
  array_finalize(&flow_field->p);
  array_finalize(&flow_field->weight);
  return 0;
}

