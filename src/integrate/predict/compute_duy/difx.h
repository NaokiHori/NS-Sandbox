#if !defined(DIFX_H)
#define DIFX_H

#include "array.h"

extern int uy_difx (
    const double c,
    const array_t * const uy,
    const double dt,
    array_t * const duy
);

#endif // DIFX_H
