#if !defined(DIFX_H)
#define DIFX_H

#include "domain.h"

extern int ux_difx(
    const domain_t * const domain,
    const double c,
    double ** const ux,
    const double dt,
    double ** const dux
);

#endif // DIFX_H
