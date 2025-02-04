#include "memory.h"
#include "array.h"

int array_init(
    const size_t nx,
    const size_t ny,
    double *** const array
) {
  double * const buffer = memory_alloc(nx * ny, sizeof(double));
  *array = memory_alloc(ny, sizeof(double *));
  for (size_t j = 0; j < ny; j++) {
    (*array)[j] = buffer + nx * j;
  }
  return 0;
}

int array_finalize(
    double *** const array
) {
  memory_free((*array)[0]);
  memory_free(*array);
  *array = NULL;
  return 0;
}

