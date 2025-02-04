#include "logger.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "./integrate.h"
#include "./integrate/decide_dt.h"
#include "./integrate/predict.h"
#include "./integrate/solve_poisson.h"
#include "./integrate/correct.h"
#include "./integrate/update_pressure.h"

int integrate(
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    double * const dt
) {
  if (0 != decide_dt(flow_field, dt)) {
    LOGGER_FAILURE("failed to find time-step size");
    goto abort;
  }
  if (0 != predict(flow_field, flow_solver, *dt)) {
    LOGGER_FAILURE("failed to predict flow field");
    goto abort;
  }
  if (0 != solve_poisson(flow_field, flow_solver, *dt)) {
    LOGGER_FAILURE("failed to solve Poisson equation to find scalar potential");
    goto abort;
  }
  if (0 != correct(flow_field, flow_solver, *dt)) {
    LOGGER_FAILURE("failed to enforce incompressibility");
    goto abort;
  }
  if (0 != update_pressure(flow_field, flow_solver)) {
    LOGGER_FAILURE("failed to update pressure field");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to update flow field");
  return 1;
}

