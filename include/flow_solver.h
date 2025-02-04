#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "domain.h" // domain_t
#include "dft/dct.h" // dct_plan_t
#include "dft/rdft.h" // rdft_plan_t
#include "tridiagonal_solver.h" // tridiagonal_solver_plan_t

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
  // y direction: tridiagonal_solver-related things
  tridiagonal_solver_plan_t * tridiagonal_solver_plan;
  double * tridiagonal_solver_l;
  double * tridiagonal_solver_c;
  double * tridiagonal_solver_u;
} poisson_solver_t;

typedef struct {
  double ** psi;
  double ** dux;
  double ** duy;
  poisson_solver_t poisson_solver;
} flow_solver_t;

extern int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalize(
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
