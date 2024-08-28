#if !defined(INTEGRATE_PREDICT_H)
#define INTEGRATE_PREDICT_H

#include "flow_field.h"
#include "flow_solver.h"

extern int predict (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

#endif // INTEGRATE_INTERNAL_H
