#if defined(TEST)

#include <math.h>
#include "array.h"
#include "./test_util.h"

static const double pi = 3.141592653589793;

double get_ux (
    const double x,
    const double y
) {
  return sin(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

double get_uy (
    const double x,
    const double y
) {
  return sin(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

double get_p (
    const double x,
    const double y
) {
  return sin(2. * pi * x / LX) * sin(2. * pi * y / LY);
}

int get_array_ux (
    array_t * const ux
) {
  for (size_t j = 1; j <= NY; j++) {
    const double y = 0.5 * (2 * j - 1) * DY;
    for (size_t i = 1; i <= NX; i++) {
      const double x = 0.5 * (2 * i - 2) * DX;
      ux[j][i] = get_ux(x, y);
    }
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    ux[j][     0] = ux[j][NX];
    ux[j][NX + 1] = ux[j][ 1];
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    ux[     0][i] = ux[NY][i];
    ux[NY + 1][i] = ux[ 1][i];
  }
  return 0;
}

int get_array_uy (
    array_t * const uy
) {
  for (size_t j = 1; j <= NY; j++) {
    const double y = 0.5 * (2 * j - 2) * DY;
    for (size_t i = 1; i <= NX; i++) {
      const double x = 0.5 * (2 * i - 1) * DX;
      uy[j][i] = get_uy(x, y);
    }
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    uy[j][     0] = uy[j][NX];
    uy[j][NX + 1] = uy[j][ 1];
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    uy[     0][i] = uy[NY][i];
    uy[NY + 1][i] = uy[ 1][i];
  }
  return 0;
}

int get_array_p (
    array_t * const p
) {
  for (size_t j = 1; j <= NY; j++) {
    const double y = 0.5 * (2 * j - 1) * DY;
    for (size_t i = 1; i <= NX; i++) {
      const double x = 0.5 * (2 * i - 1) * DX;
      p[j][i] = get_p(x, y);
    }
  }
  for (size_t j = 0; j <= NY + 1; j++) {
    p[j][     0] = p[j][NX];
    p[j][NX + 1] = p[j][ 1];
  }
  for (size_t i = 0; i <= NX + 1; i++) {
    p[     0][i] = p[NY][i];
    p[NY + 1][i] = p[ 1][i];
  }
  return 0;
}

int check_error (
    const array_t * const answer,
    const array_t * const result,
    double * const error
) {
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      error[0] += pow(answer[j][i] - result[j][i], 2.);
      error[1] = fmax(error[1], fabs(answer[j][i] - result[j][i]));
    }
  }
  error[0] = pow(error[0] / NX / NY, 1. / 2.);
  return 0;
}

#else

extern char dummy;

#endif // TEST
