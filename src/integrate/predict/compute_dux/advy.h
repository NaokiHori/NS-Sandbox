#if !defined(ADVY_H)
#define ADVY_H

#include "array.h"

extern int ux_advy (
    const array_t * const uy,
    const array_t * const ux,
    const double dt,
    array_t * const dux
);

#endif // ADVY_H
