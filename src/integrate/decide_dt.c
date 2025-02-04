#include <stdio.h>
#include <math.h> // fmin, fmax, fabs, pow
#include "logger.h"
#include "param.h"
#include "./decide_dt.h"

static const size_t ndims = 2;

static const struct {
  // Courant number
  double adv;
  // Faraday number
  double dif;
} safety_factors = {
  .adv = 0.25,
  .dif = 0.95,
};

static int decide_dt_adv(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double small = 1.e-8;
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  *dt = 1.;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = ux_imin; i <= nx; i++) {
      const double denominator = fmax(small, fabs(ux[j][i]));
      *dt = fmin(*dt, dx / denominator);
    }
  }
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double denominator = fmax(small, fabs(uy[j][i]));
      *dt = fmin(*dt, dy / denominator);
    }
  }
  *dt *= safety_factors.adv;
  return 0;
}

static int decide_dt_dif(
    const domain_t * const domain,
    double * const dt
) {
  const double dx = domain->dx;
  const double dy = domain->dy;
  *dt = Re * 0.5 / ndims * pow(fmin(dx, dy), 2.);
  *dt *= safety_factors.dif;
  return 0;
}

int decide_dt(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  double dt_adv = 0.;
  double dt_dif = 0.;
  if (0 != decide_dt_adv(domain, flow_field, &dt_adv)) {
    LOGGER_FAILURE("failed to find advective time-step constraint");
    goto abort;
  }
  if (0 != decide_dt_dif(domain, &dt_dif)) {
    LOGGER_FAILURE("failed to find diffusive time-step constraint");
    goto abort;
  }
  *dt = fmin(dt_adv, dt_dif);
  return 0;
abort:
  LOGGER_FAILURE("failed to find time-step size");
  return 1;
}

