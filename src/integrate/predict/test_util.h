#if !defined(TEST_UTIL_H)
#define TEST_UTIL_H

extern double get_ux (
    const double x,
    const double y
);

extern double get_uy (
    const double x,
    const double y
);

extern double get_p (
    const double x,
    const double y
);

extern int get_array_ux (
    array_t * const ux
);

extern int get_array_uy (
    array_t * const uy
);

extern int get_array_p (
    array_t * const p
);

extern int check_error (
    const array_t * const answer,
    const array_t * const result,
    double * const error
);

#endif // TEST_UTIL_H
