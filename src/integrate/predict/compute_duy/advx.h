#if !defined(ADVX_H)
#define ADVX_H

#include "domain.h"

extern int uy_advx(
    const domain_t * const domain,
    double ** const ux,
    double ** const uy,
    const double dt,
    double ** const duy
);

#endif // ADVX_H
