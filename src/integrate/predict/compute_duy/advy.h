#if !defined(ADVY_H)
#define ADVY_H

#include "domain.h"

extern int uy_advy(
    const domain_t * const domain,
    double ** const uy,
    const double dt,
    double ** const duy
);

#endif // ADVY_H
