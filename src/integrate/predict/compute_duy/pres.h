#if !defined(PRES_H)
#define PRES_H

#include "array.h"

extern int uy_pres (
    const array_t * const p,
    const double dt,
    array_t * const duy
);

#endif // PRES_H
