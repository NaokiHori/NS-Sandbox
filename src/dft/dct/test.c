#if defined(DCT_TEST)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "dft/dct.h"

#define MY_ASSERT(cond) \
  if (!(cond)) { \
    fprintf(stderr, "Test failed: %s (%d)\n", __func__, __LINE__); \
    return 1; \
  }

#define REPORT_SUCCESS(objective) fprintf(stderr, "Test passed: %s (%s)\n", __func__, objective);

static const double pi = 3.141592653589793238462643383;

// compute DCT type 2 naively in O(N^2)
static int naive_dct_2 (
    const size_t nitems,
    const double * xs,
    double * ys
) {
  for (size_t j = 0; j < nitems; j++) {
    double * y = ys + j;
    *y = 0.;
    for (size_t i = 0; i < nitems; i++) {
      const double phase = pi * (2. * i + 1.) * j / (2. * nitems);
      *y += 2. * xs[i] * cos(phase);
    }
  }
  return 0;
}

// compute DCT type 3 naively in O(N^2)
static int naive_dct_3 (
    const size_t nitems,
    const double * xs,
    double * ys
) {
  for (size_t j = 0; j < nitems; j++) {
    double * y = ys + j;
    *y = xs[0];
    for (size_t i = 1; i < nitems; i++) {
      const double phase = pi * (2. * j + 1.) * i / (2. * nitems);
      *y += 2. * xs[i] * cos(phase);
    }
  }
  return 0;
}

static const size_t nitems_list[] = {
  1,  2,  4,  8, 16,  32,  64, 128,  256,  512, 1024,
  3,  6, 12, 24, 48,  96, 192, 384,  768, 1536, 3072,
  5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120,
};
static const size_t repeat_for = 2;

static int test0 (
    void
) {
  for (size_t n = 0; n < sizeof(nitems_list) / sizeof(nitems_list[0]); n++) {
    const size_t nitems = nitems_list[n];
    char objective[256] = {'\0'};
    sprintf(objective, "dct2 followed by dct3 should recover original result: nitems = %4zu", nitems);
    double * buffers[] = {
      memory_alloc(nitems * repeat_for, sizeof(double)),
      memory_alloc(nitems * repeat_for, sizeof(double)),
    };
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        const double v = - 0.5 + 1. * rand() / RAND_MAX;
        buffers[0][j * nitems + i] = v;
        buffers[1][j * nitems + i] = v;
      }
    }
    dct_plan_t * plan = NULL;
    MY_ASSERT(0 == dct_init_plan(nitems, repeat_for, &plan));
    MY_ASSERT(NULL != plan);
    MY_ASSERT(0 == dct_exec_f(plan, buffers[0]));
    MY_ASSERT(0 == dct_exec_b(plan, buffers[0]));
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        // NOTE: normalise
        MY_ASSERT(fabs(buffers[1][j * nitems + i] - buffers[0][j * nitems + i] / 2. / nitems) < nitems * 1.e-14);
      }
    }
    MY_ASSERT(0 == dct_destroy_plan(&plan));
    MY_ASSERT(NULL == plan);
    memory_free(buffers[0]);
    memory_free(buffers[1]);
    REPORT_SUCCESS(objective);
  }
  return 0;
}

static int test1 (
    void
) {
  for (size_t n = 0; n < sizeof(nitems_list) / sizeof(nitems_list[0]); n++) {
    const size_t nitems = nitems_list[n];
    char objective[256] = {'\0'};
    sprintf(objective, "dct2 should yield same result as the naive one: nitems = %4zu", nitems);
    double * buffers[] = {
      memory_alloc(nitems * repeat_for, sizeof(double)),
      memory_alloc(nitems * repeat_for, sizeof(double)),
    };
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        buffers[0][j * nitems + i] = - 0.5 + 1. * rand() / RAND_MAX;
      }
    }
    dct_plan_t * plan = NULL;
    MY_ASSERT(0 == dct_init_plan(nitems, repeat_for, &plan));
    MY_ASSERT(NULL != plan);
    for (size_t j = 0; j < repeat_for; j++) {
      naive_dct_2(nitems, buffers[0] + j * nitems, buffers[1] + j * nitems);
    }
    MY_ASSERT(0 == dct_exec_f(plan, buffers[0]));
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        MY_ASSERT(fabs(buffers[1][j * nitems + i] - buffers[0][j * nitems + i]) < nitems * 1.e-13);
      }
    }
    MY_ASSERT(0 == dct_destroy_plan(&plan));
    MY_ASSERT(NULL == plan);
    memory_free(buffers[0]);
    memory_free(buffers[1]);
    REPORT_SUCCESS(objective);
  }
  return 0;
}

static int test2 (
    void
) {
  for (size_t n = 0; n < sizeof(nitems_list) / sizeof(nitems_list[0]); n++) {
    const size_t nitems = nitems_list[n];
    char objective[256] = {'\0'};
    sprintf(objective, "dct3 should yield same result as the naive one: nitems = %4zu", nitems);
    double * buffers[] = {
      memory_alloc(nitems * repeat_for, sizeof(double)),
      memory_alloc(nitems * repeat_for, sizeof(double)),
    };
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        buffers[0][j * nitems + i] = - 0.5 + 1. * rand() / RAND_MAX;
      }
    }
    dct_plan_t * plan = NULL;
    MY_ASSERT(0 == dct_init_plan(nitems, repeat_for, &plan));
    MY_ASSERT(NULL != plan);
    for (size_t j = 0; j < repeat_for; j++) {
      naive_dct_3(nitems, buffers[0] + j * nitems, buffers[1] + j * nitems);
    }
    MY_ASSERT(0 == dct_exec_b(plan, buffers[0]));
    for (size_t j = 0; j < repeat_for; j++) {
      for (size_t i = 0; i < nitems; i++) {
        MY_ASSERT(fabs(buffers[1][j * nitems + i] - buffers[0][j * nitems + i]) < nitems * 1.e-13);
      }
    }
    MY_ASSERT(0 == dct_destroy_plan(&plan));
    MY_ASSERT(NULL == plan);
    memory_free(buffers[0]);
    memory_free(buffers[1]);
    REPORT_SUCCESS(objective);
  }
  return 0;
}

int main (
    void
) {
  int retval = 0;
  retval += test0();
  retval += test1();
  retval += test2();
  return retval;
}

#else
extern char dummy;
#endif // DCT_TEST
