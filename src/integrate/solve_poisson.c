#include "logger.h"
#include "dft/rdft.h"
#include "dft/dct.h"
#include "tridiagonal_solver.h"
#include "exchange_halo.h"
#include "./solve_poisson.h"
#include "./transpose.h"

int solve_poisson(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double ** const psi = flow_solver->psi;
  double * const buf0 = poisson_solver->buf0;
  double * const buf1 = poisson_solver->buf1;
  // assign right-hand side of Poisson equation
  {
    double ** const ux = flow_field->ux;
    double ** const uy = flow_field->uy;
    const double factor = 1. / dt / poisson_solver->dft_norm;
#pragma omp parallel for
    for (size_t j = 1; j <= ny; j++) {
      for (size_t i = 1; i <= nx; i++) {
        const double dux = - ux[j    ][i    ]
                           + ux[j    ][i + 1];
        const double duy = - uy[j    ][i    ]
                           + uy[j + 1][i    ];
        const double div = (
            + 1. / dx * dux
            + 1. / dy * duy
        );
        buf0[(j - 1) * nx + (i - 1)] = factor * div;
      }
    }
  }
  // project x to wave space
  if (X_PERIODIC) {
    rdft_plan_t * const rdft_plan = poisson_solver->rdft_plan;
    if (0 != rdft_exec_f(rdft_plan, buf0)) {
      LOGGER_FAILURE("failed to perform RDFT");
      goto abort;
    }
  } else {
    dct_plan_t * const dct_plan = poisson_solver->dct_plan;
    if (0 != dct_exec_f(dct_plan, buf0)) {
      LOGGER_FAILURE("failed to perform DCT2");
      goto abort;
    }
  }
  // x-align to y-align
  if (0 != transpose(nx, ny, buf0, buf1)) {
    LOGGER_FAILURE("failed to transpose array from x-aligned to y-aligned");
    goto abort;
  }
  // solve linear systems in y
  {
    tridiagonal_solver_plan_t * const tridiagonal_solver_plan = poisson_solver->tridiagonal_solver_plan;
    const double * const tridiagonal_solver_l = poisson_solver->tridiagonal_solver_l;
    const double * const tridiagonal_solver_c = poisson_solver->tridiagonal_solver_c;
    const double * const tridiagonal_solver_u = poisson_solver->tridiagonal_solver_u;
    const double * const wavenumbers = poisson_solver->wavenumbers;
    if (0 != tridiagonal_solver_exec(tridiagonal_solver_plan, tridiagonal_solver_l, tridiagonal_solver_c, tridiagonal_solver_u, wavenumbers, buf1)) {
      LOGGER_FAILURE("failed to solve tri-diagonal matrix");
      goto abort;
    }
  }
  // y-align to x-align
  if (0 != transpose(ny, nx, buf1, buf0)) {
    LOGGER_FAILURE("failed to transpose array from y-aligned to x-aligned");
    goto abort;
  }
  // project x to physical space
  if (X_PERIODIC) {
    rdft_plan_t * const rdft_plan = poisson_solver->rdft_plan;
    if (0 != rdft_exec_b(rdft_plan, buf0)) {
      LOGGER_FAILURE("failed to perform IRDFT");
      goto abort;
    }
  } else {
    dct_plan_t * const dct_plan = poisson_solver->dct_plan;
    if (0 != dct_exec_b(dct_plan, buf0)) {
      LOGGER_FAILURE("failed to perform DCT3");
      goto abort;
    }
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      psi[j][i] = buf0[(j - 1) * nx + (i - 1)];
    }
  }
  // exchange halo
  // NOTE: since DCT assumes dpdx = 0,
  //       boundary conditions are not directly imposed
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(domain, psi)) {
      LOGGER_FAILURE("failed to exchange halo in x");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(domain, psi)) {
      LOGGER_FAILURE("failed to exchange halo in y");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

