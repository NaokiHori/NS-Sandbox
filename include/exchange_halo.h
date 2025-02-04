#if !defined(EXCHANGE_HALO_H)
#define EXCHANGE_HALO_H

#include "domain.h" // domain_t

extern int exchange_halo_x(
    const domain_t * const domain,
    double ** const array
);

extern int exchange_halo_y(
    const domain_t * const domain,
    double ** const array
);

#endif // EXCHANGE_HALO_H
