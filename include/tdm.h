#if !defined(TDM_H)
#define TDM_H

#include <stddef.h> // size_t
#include <stdbool.h> // bool

typedef struct tdm_internal_t tdm_internal_t;

typedef struct {
  // size of linear system
  size_t nitems;
  // repeat several times
  size_t repeat_for;
  bool is_periodic;
  // for internal use, opaque pointer
  tdm_internal_t * internal;
} tdm_plan_t;

extern int tdm_init_plan (
    const size_t nitems,
    const size_t repeat_for,
    const bool is_periodic,
    tdm_plan_t ** const tdm_plan
);

extern int tdm_solve (
    tdm_plan_t * const tdm_plan,
    // tri-diagonal matrix, lower, center, upper-diagonals
    const double * const l,
    const double * const c,
    const double * const u,
    // offset for center-diagonal components, can vary for each repeat
    const double * const c_offsets,
    // input and output
    double * const q
);

extern int tdm_destroy_plan (
    tdm_plan_t ** const tdm_plan
);

#endif // TDM_H
