#if !defined(MONITOR_H)
#define MONITOR_H

#include <stddef.h> // size_t
#include "flow_field.h"

extern int monitor (
    const size_t step,
    const double time,
    const double dt,
    const flow_field_t * const flow_field
);

#endif // MONITOR_H
