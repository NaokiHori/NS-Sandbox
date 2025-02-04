#if !defined(TRANSPOSE_H)
#define TRANSPOSE_H

#include <stddef.h> // size_t

extern int transpose(
    const size_t nx,
    const size_t ny,
    const double * const buf0,
    double * const buf11
);

#endif // TRANSPOSE_H
