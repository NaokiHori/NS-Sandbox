#include "array.h"
#include "domain.h"
#include "./advy.h"

int uy_advy (
    const array_t * const uy,
    const double dt,
    array_t * const duy
) {
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double uy_ym = + 0.5 * uy[j - 1][i    ]
                           + 0.5 * uy[j    ][i    ];
      const double uy_yp = + 0.5 * uy[j    ][i    ]
                           + 0.5 * uy[j + 1][i    ];
      const double duy_ym = - uy[j - 1][i    ]
                            + uy[j    ][i    ];
      const double duy_yp = - uy[j    ][i    ]
                            + uy[j + 1][i    ];
      duy[j][i] -= dt * (
          + 0.5 / DY * uy_ym * duy_ym
          + 0.5 / DY * uy_yp * duy_yp
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
  array_t * const     uy = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const result = malloc((NX + 2) * (NY + 2) * sizeof(double));
  array_t * const answer = malloc((NX + 2) * (NY + 2) * sizeof(double));
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_uy(uy);
  for (size_t j = uy_jmin; j <= NY; j++) {
    const double y = get_y(j);
    for (size_t i = 1; i <= NX; i++) {
      const double x = get_x(i);
      answer[j][i] = - get_uy(x, y) * get_duydy(x, y);
    }
  }
  uy_advy(uy, 1., result);
  double error[2] = {0., 0.};
  check_error(answer, result, error);
  printf("%6d % .15e % .15e\n", NX, error[0], error[1]);
  free(    uy);
  free(result);
  free(answer);
  return 0;
}

#endif // TEST
