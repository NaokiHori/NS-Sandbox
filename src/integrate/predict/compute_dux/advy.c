#include "array.h"
#include "domain.h"
#include "./advy.h"

int ux_advy (
    const array_t * const uy,
    const array_t * const ux,
    const double dt,
    array_t * const dux
) {
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      const double uy_ym = + 0.5 * uy[j    ][i - 1]
                           + 0.5 * uy[j    ][i    ];
      const double uy_yp = + 0.5 * uy[j + 1][i - 1]
                           + 0.5 * uy[j + 1][i    ];
      const double dux_ym = - ux[j - 1][i    ]
                            + ux[j    ][i    ];
      const double dux_yp = - ux[j    ][i    ]
                            + ux[j + 1][i    ];
      dux[j][i] -= dt * (
          + 0.5 / DY * uy_ym * dux_ym
          + 0.5 / DY * uy_yp * dux_yp
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
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_ux(ux);
  get_array_uy(uy);
  for (size_t j = 1; j <= NY; j++) {
    const double y = get_y(j);
    for (size_t i = ux_imin; i <= NX; i++) {
      const double x = get_x(i);
      answer[j][i] = - get_uy(x, y) * get_duxdy(x, y);
    }
  }
  ux_advy(uy, ux, 1., result);
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
