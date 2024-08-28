#if !defined(COMPUTE_DUY_H)
#define COMPUTE_DUY_H

#include "array.h"
#include "flow_field.h"

int compute_duy (
    const flow_field_t * const flow_field,
    const double dt,
    array_t * const duy
);

#endif // COMPUTE_DUY_H
