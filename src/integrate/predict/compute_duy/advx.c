#include "array.h"
#include "domain.h"
#include "./advx.h"

int uy_advx (
    const array_t * const ux,
    const array_t * const uy,
    const double dt,
    array_t * const duy
) {
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double ux_xm = + 0.5 * ux[j - 1][i    ]
                           + 0.5 * ux[j    ][i    ];
      const double ux_xp = + 0.5 * ux[j - 1][i + 1]
                           + 0.5 * ux[j    ][i + 1];
      const double duy_xm = - uy[j    ][i - 1]
                            + uy[j    ][i    ];
      const double duy_xp = - uy[j    ][i    ]
                            + uy[j    ][i + 1];
      duy[j][i] -= dt * (
          + 0.5 / DX * ux_xm * duy_xm
          + 0.5 / DX * ux_xp * duy_xp
      );
    }
  }
  return 0;
}

#if defined(TEST)

#include <stdio.h>
#include <stdlib.h>
#include "../test_util.h"
#include "./test_util.h"

int main (
    void
) {
  array_t * const     ux = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const     uy = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const result = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const answer = malloc((NX + 2) * (NY + 2) * sizeof(double));
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_ux(ux);
  get_array_uy(uy);
  for (size_t j = uy_jmin; j <= NY; j++) {
    const double y = get_y(j);
    for (size_t i = 1; i <= NX; i++) {
      const double x = get_x(i);
      answer[j][i] = - get_ux(x, y) * get_duydx(x, y);
    }
  }
  uy_advx(ux, uy, 1., result);
  double error[2] = {0., 0.};
  check_error(answer, result, error);
  printf("%6d % .15e % .15e\n", NX, error[0], error[1]);
  free(    ux);
  free(    uy);
  free(result);
  free(answer);
  return 0;
}

#endif // TEST
