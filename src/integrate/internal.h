#if !defined(INTEGRATE_INTERNAL_H)
#define INTEGRATE_INTERNAL_H

#include "flow_field.h"
#include "flow_solver.h"

extern int decide_dt (
    const flow_field_t * const flow_field,
    double * const dt
);

extern int solve_poisson (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

extern int correct (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

extern int update_pressure (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
);

#endif // INTEGRATE_INTERNAL_H
