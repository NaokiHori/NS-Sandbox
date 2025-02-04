#if !defined(COMPUTE_DUY_H)
#define COMPUTE_DUY_H

#include "domain.h"
#include "flow_field.h"

int compute_duy(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const double dt,
    double ** const duy
);

#endif // COMPUTE_DUY_H
