#if !defined(PRES_H)
#define PRES_H

#include "domain.h"

extern int ux_pres(
    const domain_t * const domain,
    double ** const p,
    const double dt,
    double ** const dux
);

#endif // PRES_H
