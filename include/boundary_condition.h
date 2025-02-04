#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include "domain.h" // domain_t

extern int impose_boundary_condition_ux_x(
    const domain_t * const domain,
    double ** const ux
);

extern int impose_boundary_condition_ux_y(
    const domain_t * const domain,
    double ** const ux
);

extern int impose_boundary_condition_uy_x(
    const domain_t * const domain,
    double ** const uy
);

extern int impose_boundary_condition_uy_y(
    const domain_t * const domain,
    double ** const uy
);

#endif // BOUNDARY_CONDITION_H
