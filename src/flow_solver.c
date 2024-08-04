#include <math.h>
#include "memory.h"
#include "logger.h"
#include "domain.h"
#include "flow_solver.h"
#include "dct.h"
#include "tdm.h"

static const double pi = 3.141592653589793;

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
  // Poisson solver
  double ** const buf0 = &flow_solver->poisson_solver.buf0;
  double ** const buf1 = &flow_solver->poisson_solver.buf1;
  *buf0 = memory_alloc(NX * NY, sizeof(double));
  *buf1 = memory_alloc(NX * NY, sizeof(double));
  dct_plan_t ** const dct_plan = &flow_solver->poisson_solver.dct_plan;
  tdm_plan_t ** const tdm_plan = &flow_solver->poisson_solver.tdm_plan;
  if (0 != dct_init_plan(NX, NY, dct_plan)) {
    LOGGER_FAILURE("failed to initialise DCT solver");
    goto abort;
  }
  if (0 != tdm_init_plan(NY, NX, Y_PERIODIC, tdm_plan)) {
    LOGGER_FAILURE("failed to initialise TDM solver");
    goto abort;
  }
  // DCT2 followed by DCT3 gives original array multiplied by "2 * NX"
  double ** const wave_numbers = &flow_solver->poisson_solver.wave_numbers;
  double * const dct_norm = &flow_solver->poisson_solver.dct_norm;
  *wave_numbers = memory_alloc(NX, sizeof(double));
  *dct_norm = 2. * NX;
  for (size_t i = 0; i < NX; i++) {
    (*wave_numbers)[i] = - pow(2. / DX * sin(pi * i / (2. * NX)), 2.);
  }
  double ** const tdm_l = &flow_solver->poisson_solver.tdm_l;
  double ** const tdm_c = &flow_solver->poisson_solver.tdm_c;
  double ** const tdm_u = &flow_solver->poisson_solver.tdm_u;
  *tdm_l = memory_alloc(NY, sizeof(double));
  *tdm_c = memory_alloc(NY, sizeof(double));
  *tdm_u = memory_alloc(NY, sizeof(double));
  for (size_t j = 0; j < NY; j++) {
    const double l = 1. / DY / DY;
    const double u = 1. / DY / DY;
    (*tdm_l)[j] = + 1. * l;
    (*tdm_u)[j] = + 1. * u;
    (*tdm_c)[j] = - 1. * l - 1. * u;
    // for non-periodic cases, impose Neumann boundary condition dp/dy = 0
    if (!Y_PERIODIC &&      0 == j) {
      (*tdm_c)[j] += 1. * l;
    }
    if (!Y_PERIODIC && NY - 1 == j) {
      (*tdm_c)[j] += 1. * u;
    }
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to initialise flow solver");
  return 1;
}

int flow_solver_finalise (
    flow_solver_t * const flow_solver
) {
  memory_free(flow_solver->psi);
  memory_free(flow_solver->dux);
  memory_free(flow_solver->duy);
  // Poisson solver
  memory_free(flow_solver->poisson_solver.buf0);
  memory_free(flow_solver->poisson_solver.buf1);
  dct_destroy_plan(&flow_solver->poisson_solver.dct_plan);
  tdm_destroy_plan(&flow_solver->poisson_solver.tdm_plan);
  memory_free(flow_solver->poisson_solver.wave_numbers);
  memory_free(flow_solver->poisson_solver.tdm_l);
  memory_free(flow_solver->poisson_solver.tdm_c);
  memory_free(flow_solver->poisson_solver.tdm_u);
  return 0;
}

