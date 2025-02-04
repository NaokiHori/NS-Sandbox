#if !defined(ADVX_H)
#define ADVX_H

#include "domain.h"

extern int ux_advx(
    const domain_t * const domain,
    double ** const ux,
    const double dt,
    double ** const dux
);

#endif // ADVX_H
