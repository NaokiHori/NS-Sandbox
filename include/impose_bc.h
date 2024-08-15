#if !defined(IMPOSE_BC_H)
#define IMPOSE_BC_H

#include "array.h" // array_t

extern int impose_bc_ux_x (
    array_t * const ux
);

extern int impose_bc_ux_y (
    array_t * const ux
);

extern int impose_bc_uy_x (
    array_t * const uy
);

extern int impose_bc_uy_y (
    array_t * const uy
);

#endif // IMPOSE_BC_H
