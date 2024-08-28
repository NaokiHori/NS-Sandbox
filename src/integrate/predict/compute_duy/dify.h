#if !defined(DIFY_H)
#define DIFY_H

#include "array.h"

extern int uy_dify (
    const double c,
    const array_t * const uy,
    const double dt,
    array_t * const duy
);

#endif // DIFY_H
