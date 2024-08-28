#if !defined(ADVX_H)
#define ADVX_H

#include "array.h"

extern int ux_advx (
    const array_t * const ux,
    const double dt,
    array_t * const dux
);

#endif // ADVX_H
