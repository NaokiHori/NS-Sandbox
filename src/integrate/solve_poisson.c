#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "dft/rdft.h"
#include "dft/dct.h"
#include "tdm.h"
#include "transpose.h"
#include "exchange_halo.h"
#include "./internal.h"

int solve_poisson (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  array_t * const psi = flow_solver->psi;
  double * const buf0 = poisson_solver->buf0;
  double * const buf1 = poisson_solver->buf1;
  // assign right-hand side of Poisson equation
  {
    const array_t * const ux = flow_field->ux;
    const array_t * const uy = flow_field->uy;
    const double factor = 1. / dt / poisson_solver->dft_norm;
#pragma omp parallel for
    for (size_t j = 1; j <= NY; j++) {
      for (size_t i = 1; i <= NX; i++) {
        const double dux = - ux[j    ][i    ]
                           + ux[j    ][i + 1];
        const double duy = - uy[j    ][i    ]
                           + uy[j + 1][i    ];
        const double div = (
            + 1. / DX * dux
            + 1. / DY * duy
        );
        buf0[(j - 1) * NX + (i - 1)] = factor * div;
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
  if (0 != transpose(NX, NY, buf0, buf1)) {
    LOGGER_FAILURE("failed to transpose array from x-aligned to y-aligned");
    goto abort;
  }
  // solve linear systems in y
  {
    tdm_plan_t * const tdm_plan = poisson_solver->tdm_plan;
    const double * const tdm_l = poisson_solver->tdm_l;
    const double * const tdm_c = poisson_solver->tdm_c;
    const double * const tdm_u = poisson_solver->tdm_u;
    const double * const wavenumbers = poisson_solver->wavenumbers;
    if (0 != tdm_solve(tdm_plan, tdm_l, tdm_c, tdm_u, wavenumbers, buf1)) {
      LOGGER_FAILURE("failed to solve tri-diagonal matrix");
      goto abort;
    }
  }
  // y-align to x-align
  if (0 != transpose(NY, NX, buf1, buf0)) {
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
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      psi[j][i] = buf0[(j - 1) * NX + (i - 1)];
    }
  }
  // exchange halo
  // NOTE: since DCT assumes dpdx = 0,
  //       boundary conditions are not directly imposed
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(psi)) {
      LOGGER_FAILURE("failed to exchange halo in x");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(psi)) {
      LOGGER_FAILURE("failed to exchange halo in y");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

