#if defined(TRANSPOSE_TEST)

#include <stdio.h>
#include "../transpose.h"

#define MY_ASSERT(cond) \
  if (!(cond)) { \
    fprintf(stderr, "Test failed: %s (%d)\n", __func__, __LINE__); \
    return 1; \
  }

#define REPORT_SUCCESS(objective) fprintf(stderr, "Test passed: %s (%s)\n", __func__, objective);

static int test0(
    void
) {
  const char objective[] = "compare two matrices before and after transposal";
#define nx 4
#define ny 3
  const double xs[nx * ny] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.};
  double ys[nx * ny] = {0};
  transpose(nx, ny, xs, ys);
  MY_ASSERT( 0. == ys[ 0]);
  MY_ASSERT( 4. == ys[ 1]);
  MY_ASSERT( 8. == ys[ 2]);
  MY_ASSERT( 1. == ys[ 3]);
  MY_ASSERT( 5. == ys[ 4]);
  MY_ASSERT( 9. == ys[ 5]);
  MY_ASSERT( 2. == ys[ 6]);
  MY_ASSERT( 6. == ys[ 7]);
  MY_ASSERT(10. == ys[ 8]);
  MY_ASSERT( 3. == ys[ 9]);
  MY_ASSERT( 7. == ys[10]);
  MY_ASSERT(11. == ys[11]);
  REPORT_SUCCESS(objective);
  return 0;
}

int main(
    void
) {
  int retval = 0;
  retval += test0();
  return retval;
}

#else
extern char dummy;
#endif // TRANSPOSE_TEST
