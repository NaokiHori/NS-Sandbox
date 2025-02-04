#if !defined(DECIDE_DT_H)
#define DECIDE_DT_H

#include "flow_field.h"
#include "flow_solver.h"

extern int decide_dt(
    const flow_field_t * const flow_field,
    double * const dt
);

#endif // DECIDE_DT_H
