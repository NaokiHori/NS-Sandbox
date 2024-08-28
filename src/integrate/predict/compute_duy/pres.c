#include "./pres.h"

int uy_pres (
    const array_t * const p,
    const double dt,
    array_t * const duy
) {
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      duy[j][i] -= dt / DY * (
          - p[j - 1][i    ]
          + p[j    ][i    ]
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
  array_t * const      p = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const result = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const answer = malloc((NX + 2) * (NY + 2) * sizeof(double));
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_p(p);
  for (size_t j = uy_jmin; j <= NY; j++) {
    const double y = get_y(j);
    for (size_t i = 1; i <= NX; i++) {
      const double x = get_x(i);
      answer[j][i] = - get_dpdy(x, y);
    }
  }
  uy_pres(p, 1., result);
  double error[2] = {0., 0.};
  check_error(answer, result, error);
  printf("%6d % .15e % .15e\n", NX, error[0], error[1]);
  free(     p);
  free(result);
  free(answer);
  return 0;
}

#endif // TEST
