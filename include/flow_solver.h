#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "array.h" // array_t
#include "dct.h" // dct_plan_t
#include "tdm.h" // tdm_plan_t

typedef struct {
  array_t * psi;
  array_t * dux;
  array_t * duy;
  // variables used by Poisson solver
  struct {
    // buffers to store intermediate data
    double * buf0;
    double * buf1;
    // dct-related things
    dct_plan_t * dct_plan;
    double dct_norm;
    double * wave_numbers;
    // tdm-related things
    tdm_plan_t * tdm_plan;
    double * tdm_l;
    double * tdm_c;
    double * tdm_u;
  } poisson_solver;
} flow_solver_t;

extern int flow_solver_init (
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalise (
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
