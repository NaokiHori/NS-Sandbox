#if !defined(TEST_UTIL_H)
#define TEST_UTIL_H

#include "domain.h"

extern double get_ux(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_uy(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_p(
    const domain_t * const domain,
    const double x,
    const double y
);

extern int get_array_ux(
    const domain_t * const domain,
    double ** const ux
);

extern int get_array_uy(
    const domain_t * const domain,
    double ** const uy
);

extern int get_array_p(
    const domain_t * const domain,
    double ** const p
);

extern int check_error(
    const domain_t * const domain,
    double ** const answer,
    double ** const result,
    double * const error
);

#endif // TEST_UTIL_H
