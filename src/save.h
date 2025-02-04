#if !defined(SAVE_H)
#define SAVE_H

#include <stddef.h> // size_t
#include "flow_field.h" // flow_field_t

extern int save(
    const size_t id,
    const size_t step,
    const double time,
    const flow_field_t * const flow_field
);

#endif // SAVE_H
