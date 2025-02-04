#include <math.h>
#include "memory.h"
#include "array.h"
#include "logger.h"
#include "flow_solver.h"
#include "dft/rdft.h"
#include "dft/dct.h"
#include "tridiagonal_solver.h"

static int init_x_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  double * const dft_norm = &poisson_solver->dft_norm;
  double ** const wavenumbers = &poisson_solver->wavenumbers;
  *wavenumbers = memory_alloc(nx, sizeof(double));
  if (X_PERIODIC) {
    rdft_plan_t ** const rdft_plan = &poisson_solver->rdft_plan;
    if (0 != rdft_init_plan(nx, ny, rdft_plan)) {
      LOGGER_FAILURE("failed to initialise RDFT solver");
      goto abort;
    }
    // RDFT followed by IRDFT gives original array multiplied by nx
    *dft_norm = 1. * nx;
  } else {
    dct_plan_t ** const dct_plan = &poisson_solver->dct_plan;
    if (0 != dct_init_plan(nx, ny, dct_plan)) {
      LOGGER_FAILURE("failed to initialise DCT solver");
      goto abort;
    }
    // DCT2 followed by DCT3 gives original array multiplied by 2 * nx
    *dft_norm = 2. * nx;
  }
  for (size_t i = 0; i < nx; i++) {
    const double pi = 3.1415926535897932385;
    (*wavenumbers)[i] = - pow(2. / dx * sin(pi * i / *dft_norm), 2.);
  }
  return 0;
abort:
  return 1;
}

static int init_y_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dy = domain->dy;
  tridiagonal_solver_plan_t ** const tridiagonal_solver_plan = &poisson_solver->tridiagonal_solver_plan;
  if (0 != tridiagonal_solver_init_plan(ny, nx, Y_PERIODIC, tridiagonal_solver_plan)) {
    LOGGER_FAILURE("failed to initialise tridiagonal_solver solver");
    goto abort;
  }
  double ** const tridiagonal_solver_l = &poisson_solver->tridiagonal_solver_l;
  double ** const tridiagonal_solver_c = &poisson_solver->tridiagonal_solver_c;
  double ** const tridiagonal_solver_u = &poisson_solver->tridiagonal_solver_u;
  *tridiagonal_solver_l = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_c = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_u = memory_alloc(ny, sizeof(double));
  for (size_t j = 0; j < ny; j++) {
    const double l = 1. / dy / dy;
    const double u = 1. / dy / dy;
    (*tridiagonal_solver_l)[j] = + 1. * l;
    (*tridiagonal_solver_u)[j] = + 1. * u;
    (*tridiagonal_solver_c)[j] = - 1. * l
                  - 1. * u;
    // for non-periodic cases,
    //   impose Neumann boundary condition
    //   dp/dy = 0
    if (!Y_PERIODIC &&      0 == j) {
      (*tridiagonal_solver_c)[j] += 1. * l;
    }
    if (!Y_PERIODIC && ny - 1 == j) {
      (*tridiagonal_solver_c)[j] += 1. * u;
    }
  }
  return 0;
abort:
  return 1;
}

int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  // auxiliary buffers
  array_init(nx + 2, ny + 2, &flow_solver->psi);
  array_init(nx + 2, ny + 2, &flow_solver->dux);
  array_init(nx + 2, ny + 2, &flow_solver->duy);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double ** const buf0 = &poisson_solver->buf0;
  double ** const buf1 = &poisson_solver->buf1;
  *buf0 = memory_alloc(nx * ny, sizeof(double));
  *buf1 = memory_alloc(nx * ny, sizeof(double));
  // x direction: dft-related things
  if (0 != init_x_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialise dft part of poisson solver");
    goto abort;
  }
  // y direction: tridiagonal_solver-related things
  if (0 != init_y_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialise tridiagonal_solver part of poisson solver");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to initialise flow solver");
  return 1;
}

int flow_solver_finalize(
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  array_finalize(&flow_solver->psi);
  array_finalize(&flow_solver->dux);
  array_finalize(&flow_solver->duy);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  memory_free(poisson_solver->buf0);
  memory_free(poisson_solver->buf1);
  if (X_PERIODIC) {
    rdft_destroy_plan(&poisson_solver->rdft_plan);
  } else {
    dct_destroy_plan(&poisson_solver->dct_plan);
  }
  tridiagonal_solver_destroy_plan(&poisson_solver->tridiagonal_solver_plan);
  memory_free(poisson_solver->wavenumbers);
  memory_free(poisson_solver->tridiagonal_solver_l);
  memory_free(poisson_solver->tridiagonal_solver_c);
  memory_free(poisson_solver->tridiagonal_solver_u);
  return 0;
}

