#if !defined(PRES_H)
#define PRES_H

#include "array.h"

extern int ux_pres (
    const array_t * const p,
    const double dt,
    array_t * const dux
);

#endif // PRES_H
