#if !defined(SOLVE_POISSON_H)
#define SOLVE_POISSON_H

#include "flow_field.h"
#include "flow_solver.h"

extern int solve_poisson(
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

#endif // SOLVE_POISSON_H
