#if !defined(UPDATE_PRESSURE_H)
#define UPDATE_PRESSURE_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"

extern int update_pressure(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
);

#endif // UPDATE_PRESSURE_H
