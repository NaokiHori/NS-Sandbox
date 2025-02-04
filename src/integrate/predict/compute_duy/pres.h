#if !defined(PRES_H)
#define PRES_H

#include "domain.h"

extern int uy_pres(
    const domain_t * const domain,
    double ** const p,
    const double dt,
    double ** const duy
);

#endif // PRES_H
