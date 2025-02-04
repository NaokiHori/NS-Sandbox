#include "./advy.h"

int ux_advy(
    const domain_t * const domain,
    double ** const uy,
    double ** const ux,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dy = domain->dy;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = ux_imin; i <= nx; i++) {
      const double uy_ym = + 0.5 * uy[j    ][i - 1]
                           + 0.5 * uy[j    ][i    ];
      const double uy_yp = + 0.5 * uy[j + 1][i - 1]
                           + 0.5 * uy[j + 1][i    ];
      const double dux_ym = - ux[j - 1][i    ]
                            + ux[j    ][i    ];
      const double dux_yp = - ux[j    ][i    ]
                            + ux[j + 1][i    ];
      dux[j][i] -= dt * (
          + 0.5 / dy * uy_ym * dux_ym
          + 0.5 / dy * uy_yp * dux_yp
      );
    }
  }
  return 0;
}

#if defined(TEST)

#include <stdio.h> // printf
#include <stdlib.h> // strtol
#include "array.h"
#include "domain.h"
#include "../test_util.h"
#include "./test_util.h"

int main(
    int argc,
    char * argv[]
) {
  if (2 != argc) {
    printf("invalid number of arguments: %d, expected 2\n", argc);
    return 1;
  }
  const double length = 1.;
  const size_t nx = strtol(argv[1], NULL, 10);
  const size_t ny = strtol(argv[1], NULL, 10);
  const domain_t domain = {
    .lx = length,
    .ly = length,
    .nx = nx,
    .ny = ny,
    .dx = length / nx,
    .dy = length / ny,
  };
  double ** ux = NULL;
  double ** uy = NULL;
  double ** result = NULL;
  double ** answer = NULL;
  array_init(nx + 2, ny + 2, &ux);
  array_init(nx + 2, ny + 2, &uy);
  array_init(nx + 2, ny + 2, &result);
  array_init(nx + 2, ny + 2, &answer);
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = ux_imin; i <= nx; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_ux(&domain, ux);
  get_array_uy(&domain, uy);
  for (size_t j = 1; j <= ny; j++) {
    const double y = get_y(&domain, j);
    for (size_t i = ux_imin; i <= nx; i++) {
      const double x = get_x(&domain, i);
      answer[j][i] = - get_uy(&domain, x, y) * get_duxdy(&domain, x, y);
    }
  }
  ux_advy(&domain, uy, ux, 1., result);
  double error[2] = {0., 0.};
  check_error(&domain, answer, result, error);
  printf("%6zu % .15e % .15e\n", nx, error[0], error[1]);
  array_finalize(&ux);
  array_finalize(&uy);
  array_finalize(&result);
  array_finalize(&answer);
  return 0;
}

#endif // TEST
