#if !defined(DIFY_H)
#define DIFY_H

#include "domain.h"

extern int uy_dify(
    const domain_t * const domain,
    const double c,
    double ** const uy,
    const double dt,
    double ** const duy
);

#endif // DIFY_H
