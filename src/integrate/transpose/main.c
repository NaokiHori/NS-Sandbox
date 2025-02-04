#include "../transpose.h"

int transpose(
    const size_t nx,
    const size_t ny,
    const double * const buf0,
    double * const buf1
) {
#pragma omp parallel for
  for (size_t j = 0; j < ny; j++) {
    for (size_t i = 0; i < nx; i++) {
      buf1[i * ny + j] = buf0[j * nx + i];
    }
  }
  return 0;
}

