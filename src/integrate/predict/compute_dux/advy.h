#if !defined(ADVY_H)
#define ADVY_H

#include "domain.h"

extern int ux_advy(
    const domain_t * const domain,
    double ** const uy,
    double ** const ux,
    const double dt,
    double ** const dux
);

#endif // ADVY_H
