#if !defined(ARRAY_H)
#define ARRAY_H

#include <stddef.h>

extern int array_init(
    const size_t nx,
    const size_t ny,
    double *** const array
);

extern int array_finalize(
    double *** const array
);

#endif // ARRAY_H
