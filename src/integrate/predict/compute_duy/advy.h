#if !defined(ADVY_H)
#define ADVY_H

#include "array.h"

extern int uy_advy (
    const array_t * const uy,
    const double dt,
    array_t * const duy
);

#endif // ADVY_H
