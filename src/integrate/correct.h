#if !defined(CORRECT_H)
#define CORRECT_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"

extern int correct(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt
);

#endif // CORRECT_H
