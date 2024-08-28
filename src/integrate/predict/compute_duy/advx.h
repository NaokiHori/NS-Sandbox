#if !defined(ADVX_H)
#define ADVX_H

#include "array.h"

extern int uy_advx (
    const array_t * const ux,
    const array_t * const uy,
    const double dt,
    array_t * const duy
);

#endif // ADVX_H
