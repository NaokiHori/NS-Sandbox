#if !defined(MONITOR_H)
#define MONITOR_H

#include <stddef.h> // size_t
#include "domain.h" // domain_t
#include "flow_field.h" // flow_field_t

extern int monitor(
    const size_t step,
    const double time,
    const double dt,
    const domain_t * const domain,
    const flow_field_t * const flow_field
);

#endif // MONITOR_H
