#if !defined(DIFY_H)
#define DIFY_H

#include "domain.h"

extern int ux_dify(
    const domain_t * const domain,
    const double c,
    double ** const ux,
    const double dt,
    double ** const dux
);

#endif // DIFY_H
