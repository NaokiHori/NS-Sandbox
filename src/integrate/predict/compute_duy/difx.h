#if !defined(DIFX_H)
#define DIFX_H

#include "domain.h"

extern int uy_difx(
    const domain_t * const domain,
    const double c,
    double ** const uy,
    const double dt,
    double ** const duy
);

#endif // DIFX_H
