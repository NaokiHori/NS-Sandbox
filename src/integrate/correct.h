#if !defined(CORRECT_H)
#define CORRECT_H

#include "flow_field.h"
#include "flow_solver.h"

extern int correct(
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

#endif // CORRECT_H
