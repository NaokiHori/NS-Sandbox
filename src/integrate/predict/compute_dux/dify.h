#if !defined(DIFY_H)
#define DIFY_H

#include "array.h"

extern int ux_dify (
    const double c,
    const array_t * const ux,
    const double dt,
    array_t * const dux
);

#endif // DIFY_H
