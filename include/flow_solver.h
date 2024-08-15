#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "array.h" // array_t
#include "dft/dct.h" // dct_plan_t
#include "dft/rdft.h" // rdft_plan_t
#include "tdm.h" // tdm_plan_t

// variables used to solve Poisson equations
typedef struct {
  // buffers to store intermediate data
  double * buf0;
  double * buf1;
  // x direction: dft-related things
  // NOTE: rdft_plan is used when X_PERIODIC; otherwise dct_plan is used
  rdft_plan_t * rdft_plan;
  dct_plan_t * dct_plan;
  double dft_norm;
  double * wavenumbers;
  // y direction: tdm-related things
  tdm_plan_t * tdm_plan;
  double * tdm_l;
  double * tdm_c;
  double * tdm_u;
} poisson_solver_t;

typedef struct {
  array_t * psi;
  array_t * dux;
  array_t * duy;
  poisson_solver_t poisson_solver;
} flow_solver_t;

extern int flow_solver_init (
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalise (
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
