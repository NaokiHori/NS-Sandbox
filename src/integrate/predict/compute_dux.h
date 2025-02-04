#if !defined(COMPUTE_DUX_H)
#define COMPUTE_DUX_H

#include "domain.h"
#include "flow_field.h"

int compute_dux(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const double dt,
    double ** const dux
);

#endif // COMPUTE_DUX_H
