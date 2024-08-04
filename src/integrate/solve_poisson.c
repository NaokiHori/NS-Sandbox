#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "dct.h"
#include "tdm.h"
#include "transpose.h"
#include "impose_bc.h"
#include "exchange_halo.h"
#include "./internal.h"

int solve_poisson (
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  if (X_PERIODIC) {
    LOGGER_FAILURE("poisson solver is not ready for periodic boundary in x");
    goto abort;
  }
  array_t * const psi = flow_solver->psi;
  double * const buf0 = flow_solver->poisson_solver.buf0;
  double * const buf1 = flow_solver->poisson_solver.buf1;
  // assign right-hand side of Poisson equation
  {
    const array_t * const ux = flow_field->ux;
    const array_t * const uy = flow_field->uy;
    const double factor = 1. / dt / flow_solver->poisson_solver.dct_norm;
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
  {
    dct_plan_t * const dct_plan = flow_solver->poisson_solver.dct_plan;
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
    tdm_plan_t * const tdm_plan = flow_solver->poisson_solver.tdm_plan;
    const double * const tdm_l = flow_solver->poisson_solver.tdm_l;
    const double * const tdm_c = flow_solver->poisson_solver.tdm_c;
    const double * const tdm_u = flow_solver->poisson_solver.tdm_u;
    const double * const wave_numbers = flow_solver->poisson_solver.wave_numbers;
    if (0 != tdm_solve(tdm_plan, tdm_l, tdm_c, tdm_u, wave_numbers, buf1)) {
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
  {
    dct_plan_t * const dct_plan = flow_solver->poisson_solver.dct_plan;
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
  // exchange halo / impose bc
  if (X_PERIODIC) {
    if (0 != exchange_halo_x(psi)) {
      LOGGER_FAILURE("failed to exchange halo in x");
      goto abort;
    }
  } else {
    if (0 != impose_bc_p_x(psi)) {
      LOGGER_FAILURE("failed to impose boundary condition in x");
      goto abort;
    }
  }
  if (Y_PERIODIC) {
    if (0 != exchange_halo_y(psi)) {
      LOGGER_FAILURE("failed to exchange halo in y");
      goto abort;
    }
  } else {
    if (0 != impose_bc_p_y(psi)) {
      LOGGER_FAILURE("failed to impose boundary condition in y");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

