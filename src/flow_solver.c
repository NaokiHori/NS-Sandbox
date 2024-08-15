#include <math.h>
#include "memory.h"
#include "logger.h"
#include "domain.h"
#include "flow_solver.h"
#include "dft/rdft.h"
#include "dft/dct.h"
#include "tdm.h"

static int init_x_solver (
    poisson_solver_t * const poisson_solver
) {
  double * const dft_norm = &poisson_solver->dft_norm;
  double ** const wavenumbers = &poisson_solver->wavenumbers;
  *wavenumbers = memory_alloc(NX, sizeof(double));
  if (X_PERIODIC) {
    rdft_plan_t ** const rdft_plan = &poisson_solver->rdft_plan;
    if (0 != rdft_init_plan(NX, NY, rdft_plan)) {
      LOGGER_FAILURE("failed to initialise RDFT solver");
      goto abort;
    }
    // RDFT followed by IRDFT gives original array multiplied by NX
    *dft_norm = 1. * NX;
  } else {
    dct_plan_t ** const dct_plan = &poisson_solver->dct_plan;
    if (0 != dct_init_plan(NX, NY, dct_plan)) {
      LOGGER_FAILURE("failed to initialise DCT solver");
      goto abort;
    }
    // DCT2 followed by DCT3 gives original array multiplied by 2 * NX
    *dft_norm = 2. * NX;
  }
  for (size_t i = 0; i < NX; i++) {
    const double pi = 3.1415926535897932385;
    (*wavenumbers)[i] = - pow(2. / DX * sin(pi * i / *dft_norm), 2.);
  }
  return 0;
abort:
  return 1;
}

static int init_y_solver (
    poisson_solver_t * const poisson_solver
) {
  tdm_plan_t ** const tdm_plan = &poisson_solver->tdm_plan;
  if (0 != tdm_init_plan(NY, NX, Y_PERIODIC, tdm_plan)) {
    LOGGER_FAILURE("failed to initialise TDM solver");
    goto abort;
  }
  double ** const tdm_l = &poisson_solver->tdm_l;
  double ** const tdm_c = &poisson_solver->tdm_c;
  double ** const tdm_u = &poisson_solver->tdm_u;
  *tdm_l = memory_alloc(NY, sizeof(double));
  *tdm_c = memory_alloc(NY, sizeof(double));
  *tdm_u = memory_alloc(NY, sizeof(double));
  for (size_t j = 0; j < NY; j++) {
    const double l = 1. / DY / DY;
    const double u = 1. / DY / DY;
    (*tdm_l)[j] = + 1. * l;
    (*tdm_u)[j] = + 1. * u;
    (*tdm_c)[j] = - 1. * l
                  - 1. * u;
    // for non-periodic cases,
    //   impose Neumann boundary condition
    //   dp/dy = 0
    if (!Y_PERIODIC &&      0 == j) {
      (*tdm_c)[j] += 1. * l;
    }
    if (!Y_PERIODIC && NY - 1 == j) {
      (*tdm_c)[j] += 1. * u;
    }
  }
  return 0;
abort:
  return 1;
}

int flow_solver_init (
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  array_t ** const psi = &flow_solver->psi;
  array_t ** const dux = &flow_solver->dux;
  array_t ** const duy = &flow_solver->duy;
  *psi = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  *dux = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  *duy = memory_alloc((NX + 2) * (NY + 2), sizeof(double));
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double ** const buf0 = &poisson_solver->buf0;
  double ** const buf1 = &poisson_solver->buf1;
  *buf0 = memory_alloc(NX * NY, sizeof(double));
  *buf1 = memory_alloc(NX * NY, sizeof(double));
  // x direction: dft-related things
  if (0 != init_x_solver(poisson_solver)) {
    LOGGER_FAILURE("failed to initialise dft part of poisson solver");
    goto abort;
  }
  // y direction: tdm-related things
  if (0 != init_y_solver(poisson_solver)) {
    LOGGER_FAILURE("failed to initialise tdm part of poisson solver");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to initialise flow solver");
  return 1;
}

int flow_solver_finalise (
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  memory_free(flow_solver->psi);
  memory_free(flow_solver->dux);
  memory_free(flow_solver->duy);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  memory_free(poisson_solver->buf0);
  memory_free(poisson_solver->buf1);
  if (X_PERIODIC) {
    rdft_destroy_plan(&poisson_solver->rdft_plan);
  } else {
    dct_destroy_plan(&poisson_solver->dct_plan);
  }
  tdm_destroy_plan(&poisson_solver->tdm_plan);
  memory_free(poisson_solver->wavenumbers);
  memory_free(poisson_solver->tdm_l);
  memory_free(poisson_solver->tdm_c);
  memory_free(poisson_solver->tdm_u);
  return 0;
}

