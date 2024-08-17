#if !defined(EXCHANGE_HALO_H)
#define EXCHANGE_HALO_H

#include "array.h" // array_t

extern int exchange_halo_x (
    array_t * const array
);

extern int exchange_halo_y (
    array_t * const array
);

#endif // EXCHANGE_HALO_H
