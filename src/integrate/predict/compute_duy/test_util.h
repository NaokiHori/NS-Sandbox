#if !defined(COMPUTE_DUY_TEST_UTIL_H)
#define COMPUTE_DUY_TEST_UTIL_H

#include <stddef.h> // size_t
#include "domain.h"

extern double get_x(
    const domain_t * const domain,
    const size_t i
);

extern double get_y(
    const domain_t * const domain,
    const size_t j
);

extern double get_duydx(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_duydy(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_d2uydx2(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_d2uydy2(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_dpdy(
    const domain_t * const domain,
    const double x,
    const double y
);

#endif // COMPUTE_DUY_TEST_UTIL_H
