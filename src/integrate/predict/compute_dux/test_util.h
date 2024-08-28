#if !defined(COMPUTE_DUX_TEST_UTIL_H)
#define COMPUTE_DUX_TEST_UTIL_H

#include <stddef.h> // size_t

extern double get_x (
    const size_t i
);

extern double get_y (
    const size_t j
);

extern double get_duxdx (
    const double x,
    const double y
);

extern double get_duxdy (
    const double x,
    const double y
);

extern double get_d2uxdx2 (
    const double x,
    const double y
);

extern double get_d2uxdy2 (
    const double x,
    const double y
);

extern double get_dpdx (
    const double x,
    const double y
);

#endif // COMPUTE_DUX_TEST_UTIL_H
