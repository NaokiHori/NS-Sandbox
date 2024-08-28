#if !defined(COMPUTE_DUX_H)
#define COMPUTE_DUX_H

#include "array.h"
#include "flow_field.h"

int compute_dux (
    const flow_field_t * const flow_field,
    const double dt,
    array_t * const dux
);

#endif // COMPUTE_DUX_H
