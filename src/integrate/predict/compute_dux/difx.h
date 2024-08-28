#if !defined(DIFX_H)
#define DIFX_H

#include "array.h"

extern int ux_difx (
    const double c,
    const array_t * const ux,
    const double dt,
    array_t * const dux
);

#endif // DIFX_H
