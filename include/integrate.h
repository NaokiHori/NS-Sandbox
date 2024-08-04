#if !defined(INTEGRATE_H)
#define INTEGRATE_H

#include "flow_field.h" // flow_field_t
#include "flow_solver.h" // flow_solver_t

extern int integrate (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    double * const dt
);

#endif // INTEGRATE_H
