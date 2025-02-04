#include <math.h>
#include "memory.h"
#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "boundary_condition.h"
#include "exchange_halo.h"

static int allocate (
    flow_field_t * const flow_field
) {
  array_t ** const     ux = &flow_field->    ux;
  array_t ** const     uy = &flow_field->    uy;
  array_t ** const      p = &flow_field->     p;
  array_t ** const weight = &flow_field->weight;
  *    ux = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  *    uy = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  *     p = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  *weight = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  return 0;
}

static int init_ux (
    array_t * const ux
) {
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      ux[j][i] = 0.;
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(ux)) {
      LOGGER_FAILURE("failed to exchange halo in x (ux)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_ux_x(ux)) {
      LOGGER_FAILURE("failed to impose boundary condition in x (ux)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(ux)) {
      LOGGER_FAILURE("failed to exchange halo in y (ux)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_ux_y(ux)) {
      LOGGER_FAILURE("failed to impose boundary condition in y (ux)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_uy (
    array_t * const uy
) {
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      uy[j][i] = -1.;
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(uy)) {
      LOGGER_FAILURE("failed to exchange halo in x (uy)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_uy_x(uy)) {
      LOGGER_FAILURE("failed to impose boundary condition in x (uy)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(uy)) {
      LOGGER_FAILURE("failed to exchange halo in y (uy)");
      goto abort;
    }
  } else {
    if (0 != impose_boundary_condition_uy_y(uy)) {
      LOGGER_FAILURE("failed to impose boundary condition in y (uy)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_p (
    array_t * const p
) {
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      p[j][i] = 0.;
    }
  }
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(p)) {
      LOGGER_FAILURE("failed to exchange halo in x (p)");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(p)) {
      LOGGER_FAILURE("failed to exchange halo in y (p)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_weight (
    array_t * const weight
) {
  const double center[2] = {
    0.501 * LX,
    5. * LY / 6.,
  };
  const double radius = LX / 16.;
  for (size_t j = 0; j <= NY + 1; j++) {
    const double y = 0.5 * (2 * j - 1) * DY;
    for (size_t i = 0; i <= NX + 1; i++) {
      const double x = 0.5 * (2 * i - 1) * DX;
      const double d = sqrt(
          + pow(x - center[0], 2.)
          + pow(y - center[1], 2.)
      );
      weight[j][i] = 0.5 * (1. + tanh(NY * (d - radius)));
    }
  }
  return 0;
}

int flow_field_init (
    flow_field_t * const flow_field
) {
  if (0 != allocate(flow_field)) {
    LOGGER_FAILURE("failed to allocate flow field");
    goto abort;
  }
  if (0 != init_ux(flow_field->ux)) {
    LOGGER_FAILURE("failed to initialise ux");
    goto abort;
  }
  if (0 != init_uy(flow_field->uy)) {
    LOGGER_FAILURE("failed to initialise uy");
    goto abort;
  }
  if (0 != init_p(flow_field->p)) {
    LOGGER_FAILURE("failed to initialise p");
    goto abort;
  }
  if (0 != init_weight(flow_field->weight)) {
    LOGGER_FAILURE("failed to initialise weight");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int flow_field_finalise (
    flow_field_t * const flow_field
) {
  memory_free(flow_field->    ux);
  memory_free(flow_field->    uy);
  memory_free(flow_field->     p);
  memory_free(flow_field->weight);
  return 0;
}

