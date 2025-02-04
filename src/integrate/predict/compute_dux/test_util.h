#if !defined(COMPUTE_DUX_TEST_UTIL_H)
#define COMPUTE_DUX_TEST_UTIL_H

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

extern double get_duxdx(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_duxdy(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_d2uxdx2(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_d2uxdy2(
    const domain_t * const domain,
    const double x,
    const double y
);

extern double get_dpdx(
    const domain_t * const domain,
    const double x,
    const double y
);

#endif // COMPUTE_DUX_TEST_UTIL_H
