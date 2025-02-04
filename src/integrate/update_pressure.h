#if !defined(UPDATE_PRESSURE_H)
#define UPDATE_PRESSURE_H

#include "flow_field.h"
#include "flow_solver.h"

extern int update_pressure(
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
);

#endif // UPDATE_PRESSURE_H
