#if defined(TEST)

#include <math.h>
#include "./test_util.h"

static const double pi = 3.141592653589793;

double get_ux(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return sin(2. * pi * x / lx) * sin(2. * pi * y / ly);
}

double get_uy(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return sin(2. * pi * x / lx) * sin(2. * pi * y / ly);
}

double get_p(
    const domain_t * const domain,
    const double x,
    const double y
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  return sin(2. * pi * x / lx) * sin(2. * pi * y / ly);
}

int get_array_ux(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  for (size_t j = 1; j <= ny; j++) {
    const double y = 0.5 * (2 * j - 1) * dy;
    for (size_t i = 1; i <= nx; i++) {
      const double x = 0.5 * (2 * i - 2) * dx;
      ux[j][i] = get_ux(domain, x, y);
    }
  }
  for (size_t j = 0; j <= ny + 1; j++) {
    ux[j][     0] = ux[j][nx];
    ux[j][nx + 1] = ux[j][ 1];
  }
  for (size_t i = 0; i <= nx + 1; i++) {
    ux[     0][i] = ux[ny][i];
    ux[ny + 1][i] = ux[ 1][i];
  }
  return 0;
}

int get_array_uy(
    const domain_t * const domain,
    double ** const uy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  for (size_t j = 1; j <= ny; j++) {
    const double y = 0.5 * (2 * j - 2) * dy;
    for (size_t i = 1; i <= nx; i++) {
      const double x = 0.5 * (2 * i - 1) * dx;
      uy[j][i] = get_uy(domain, x, y);
    }
  }
  for (size_t j = 0; j <= ny + 1; j++) {
    uy[j][     0] = uy[j][nx];
    uy[j][nx + 1] = uy[j][ 1];
  }
  for (size_t i = 0; i <= nx + 1; i++) {
    uy[     0][i] = uy[ny][i];
    uy[ny + 1][i] = uy[ 1][i];
  }
  return 0;
}

int get_array_p(
    const domain_t * const domain,
    double ** const p
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  for (size_t j = 1; j <= ny; j++) {
    const double y = 0.5 * (2 * j - 1) * dy;
    for (size_t i = 1; i <= nx; i++) {
      const double x = 0.5 * (2 * i - 1) * dx;
      p[j][i] = get_p(domain, x, y);
    }
  }
  for (size_t j = 0; j <= ny + 1; j++) {
    p[j][     0] = p[j][nx];
    p[j][nx + 1] = p[j][ 1];
  }
  for (size_t i = 0; i <= nx + 1; i++) {
    p[     0][i] = p[ny][i];
    p[ny + 1][i] = p[ 1][i];
  }
  return 0;
}

int check_error(
    const domain_t * const domain,
    double ** const answer,
    double ** const result,
    double * const error
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      error[0] += pow(answer[j][i] - result[j][i], 2.);
      error[1] = fmax(error[1], fabs(answer[j][i] - result[j][i]));
    }
  }
  error[0] = pow(error[0] / nx / ny, 1. / 2.);
  return 0;
}

#else

extern char dummy;

#endif // TEST
