#include "./dify.h"

int ux_dify (
    const double c,
    const array_t * const ux,
    const double dt,
    array_t * const dux
) {
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      dux[j][i] += dt * c / DY / DY * (
          + 1. * ux[j - 1][i    ]
          - 2. * ux[j    ][i    ]
          + 1. * ux[j + 1][i    ]
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
  array_t * const result = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const answer = malloc((NX + 2) * (NY + 2) * sizeof(double));
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_ux(ux);
  for (size_t j = 1; j <= NY; j++) {
    const double y = get_y(j);
    for (size_t i = ux_imin; i <= NX; i++) {
      const double x = get_x(i);
      answer[j][i] = get_d2uxdy2(x, y);
    }
  }
  ux_dify(1., ux, 1., result);
  double error[2] = {0., 0.};
  check_error(answer, result, error);
  printf("%6d % .15e % .15e\n", NX, error[0], error[1]);
  free(    ux);
  free(result);
  free(answer);
  return 0;
}

#endif // TEST
