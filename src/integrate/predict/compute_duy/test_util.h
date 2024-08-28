#if !defined(COMPUTE_DUY_TEST_UTIL_H)
#define COMPUTE_DUY_TEST_UTIL_H

#include <stddef.h> // size_t

extern double get_x (
    const size_t i
);

extern double get_y (
    const size_t j
);

extern double get_duydx (
    const double x,
    const double y
);

extern double get_duydy (
    const double x,
    const double y
);

extern double get_d2uydx2 (
    const double x,
    const double y
);

extern double get_d2uydy2 (
    const double x,
    const double y
);

extern double get_dpdy (
    const double x,
    const double y
);

#endif // COMPUTE_DUY_TEST_UTIL_H
