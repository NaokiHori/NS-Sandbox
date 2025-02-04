#include "./pres.h"

int uy_pres(
    const domain_t * const domain,
    double ** const p,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dy = domain->dy;
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      duy[j][i] -= dt / dy * (
          - p[j - 1][i    ]
          + p[j    ][i    ]
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
  double ** p = NULL;
  double ** result = NULL;
  double ** answer = NULL;
  array_init(nx + 2, ny + 2, &p);
  array_init(nx + 2, ny + 2, &result);
  array_init(nx + 2, ny + 2, &answer);
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      result[j][i] = 0.;
    }
  }
  get_array_p(&domain, p);
  for (size_t j = uy_jmin; j <= ny; j++) {
    const double y = get_y(&domain, j);
    for (size_t i = 1; i <= nx; i++) {
      const double x = get_x(&domain, i);
      answer[j][i] = - get_dpdy(&domain, x, y);
    }
  }
  uy_pres(&domain, p, 1., result);
  double error[2] = {0., 0.};
  check_error(&domain, answer, result, error);
  printf("%6zu % .15e % .15e\n", nx, error[0], error[1]);
  array_finalize(&p);
  array_finalize(&result);
  array_finalize(&answer);
  return 0;
}

#endif // TEST
