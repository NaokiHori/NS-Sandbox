#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include "array.h" // array_t

extern int impose_boundary_condition_ux_x (
    array_t * const ux
);

extern int impose_boundary_condition_ux_y (
    array_t * const ux
);

extern int impose_boundary_condition_uy_x (
    array_t * const uy
);

extern int impose_boundary_condition_uy_y (
    array_t * const uy
);

#endif // BOUNDARY_CONDITION_H
