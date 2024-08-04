#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dct.h"

// discrete cosine transforms of type 2 and 3, Lee 1984

static const double pi = 3.141592653589793;
static const double sqrt2h = 0.7071067811865475;
static const double sqrt3h = 0.8660254037844386;

struct dct_plan_t {
  // size of the input / output signals
  size_t nitems;
  // repeat the same DCTs
  size_t repeat_for;
  // trigonometric table
  //   1 / (2 cos( (pi i) / (2 N) ))
  //   where i = 0, 1, ..., N - 1
  double * table;
  // internal buffer
  double * buf;
};

static void * memory_alloc (
    const size_t size
) {
  void * const ptr = malloc(size);
  if (NULL == ptr) {
    fprintf(stderr, "[dct FATAL] failed to allocate %zu bytes\n", size);
    return NULL;
  }
  return ptr;
}

static void memory_free (
    void * const ptr
) {
  free(ptr);
}

static int dct2 (
    const size_t nitems,
    const size_t inv,
    const double * const restrict table,
    double * const restrict xs,
    double * const restrict ys
) {
  if (1 == nitems) {
  } else if (2 == nitems) {
		const double v0 = xs[0];
		const double v1 = xs[1];
		xs[0] = 1.     * v0 + 1.     * v1;
		xs[1] = sqrt2h * v0 - sqrt2h * v1;
  } else if (3 == nitems) {
    const double v0 = xs[0];
    const double v1 = xs[1];
    const double v2 = xs[2];
    xs[0] = 1.     * v0 + 1. * v1 + 1.     * v2;
    xs[1] = sqrt3h * v0           - sqrt3h * v2;
    xs[2] = 0.5    * v0 - 1. * v1 + 0.5    * v2;
  } else if (0 == nitems % 2) {
    const size_t nhalfs = nitems / 2;
    for (size_t i = 0; i < nhalfs; i++) {
      const double c = table[(2 * i + 1) * inv];
      const double v0 = xs[             i];
      const double v1 = xs[nitems - 1 - i];
      ys[i         ] = 1. * (v0 + v1);
      ys[i + nhalfs] = c  * (v0 - v1);
    }
    dct2(nhalfs, inv * 2, table, ys +      0, xs);
    dct2(nhalfs, inv * 2, table, ys + nhalfs, xs);
    for (size_t i = 0; i < nhalfs - 1; i++) {
      xs[i * 2 + 0] = ys[         i    ];
      xs[i * 2 + 1] = ys[nhalfs + i    ]
                    + ys[nhalfs + i + 1];
    }
    xs[nitems - 2] = ys[nhalfs - 1];
    xs[nitems - 1] = ys[nitems - 1];
  } else {
    // fallback to N^2 DCT2
    for (size_t j = 0; j < nitems; j++) {
      double * y = ys + j;
      *y = 0.;
      for (size_t i = 0; i < nitems; i++) {
        const double phase = pi * (2. * i + 1.) * j / (2. * nitems);
        *y += xs[i] * cos(phase);
      }
    }
    for (size_t i = 0; i < nitems; i++) {
      xs[i] = ys[i];
    }
  }
  return 0;
}

static int dct3 (
    const size_t nitems,
    const size_t inv,
    const double * const restrict table,
    double * const restrict xs,
    double * const restrict ys
) {
  if (1 == nitems) {
  } else if (2 == nitems) {
    const double v0 = xs[0];
    const double v1 = xs[1];
    xs[0] = v0 + sqrt2h * v1;
    xs[1] = v0 - sqrt2h * v1;
  } else if (3 == nitems) {
    const double v0 = xs[0];
    const double v1 = xs[1];
    const double v2 = xs[2];
    xs[0] = v0 + sqrt3h * v1 + 0.5 * v2;
    xs[1] = v0               -       v2;
    xs[2] = v0 - sqrt3h * v1 + 0.5 * v2;
  } else if (0 == nitems % 2) {
    const size_t nhalfs = nitems / 2;
    ys[     0] = xs[0];
    ys[nhalfs] = xs[1];
    for (size_t i = 1; i < nhalfs; i++) {
      ys[         i] = xs[i * 2    ];
      ys[nhalfs + i] = xs[i * 2 - 1]
                     + xs[i * 2 + 1];
    }
    dct3(nhalfs, inv * 2, table, ys         , xs);
    dct3(nhalfs, inv * 2, table, ys + nhalfs, xs);
    for (size_t i = 0; i < nhalfs; i++) {
      const double c = table[(2 * i + 1) * inv];
      const double v0 = 1. * ys[         i];
      const double v1 = c  * ys[nhalfs + i];
      xs[             i] = v0 + v1;
      xs[nitems - 1 - i] = v0 - v1;
    }
  } else {
    // fallback to N^2 DCT3
    for (size_t j = 0; j < nitems; j++) {
      double * y = ys + j;
      *y = xs[0];
      for (size_t i = 1; i < nitems; i++) {
        const double phase = pi * (2. * j + 1.) * i / (2. * nitems);
        *y += xs[i] * cos(phase);
      }
    }
    for (size_t i = 0; i < nitems; i++) {
      xs[i] = ys[i];
    }
  }
  return 0;
}

// allocate, initalise, and pack
int dct_init_plan (
    const size_t nitems,
    const size_t repeat_for,
    dct_plan_t ** const plan
) {
  *plan = memory_alloc(1 * sizeof(dct_plan_t));
  if (NULL == *plan) {
    return 1;
  }
  (*plan)->nitems = nitems;
  (*plan)->repeat_for = repeat_for;
  // trigonometric table
  double ** table = &(*plan)->table;
  *table = memory_alloc(nitems * sizeof(double));
  if (NULL == *table) {
    return 1;
  }
  for (size_t i = 0; i < nitems; i++) {
    const double phase = (pi * i) / (2. * nitems);
    (*table)[i] = 0.5 / cos(phase);
  }
  // internal buffer
  double ** buf = &(*plan)->buf;
  *buf = memory_alloc(nitems * repeat_for * sizeof(double));
  return 0;
}

int dct_destroy_plan (
    dct_plan_t ** const plan
) {
  if (NULL == *plan) {
    fprintf(stderr, "the plan is NULL\n");
    return 1;
  }
  memory_free((*plan)->table);
  memory_free((*plan)->buf);
  memory_free(*plan);
  *plan = NULL;
  return 0;
}

int dct_exec_f (
    dct_plan_t * const plan,
    double * restrict const xs
){
  if (NULL == plan) {
    fprintf(stderr, "the plan is NULL\n");
    return 1;
  }
  const size_t nitems = plan->nitems;
  const size_t repeat_for = plan->repeat_for;
  const double * const table = plan->table;
  double * const ys = plan->buf;
#pragma omp parallel for
  for (size_t j = 0; j < repeat_for; j++) {
    dct2(nitems, 1, table, xs + j * nitems, ys + j * nitems);
    for (size_t i = 0; i < nitems; i++) {
      xs[j * nitems + i] *= 2.;
    }
  }
  return 0;
}

int dct_exec_b (
    dct_plan_t * const plan,
    double * restrict const xs
) {
  if (NULL == plan) {
    fprintf(stderr, "the plan is NULL\n");
    return 1;
  }
  const size_t nitems = plan->nitems;
  const size_t repeat_for = plan->repeat_for;
  const double * const table = plan->table;
  double * const ys = plan->buf;
#pragma omp parallel for
  for (size_t j = 0; j < repeat_for; j++) {
    xs[j * nitems + 0] *= 0.5;
    dct3(nitems, 1, table, xs + j * nitems, ys + j * nitems);
    for (size_t i = 0; i < nitems; i++) {
      xs[j * nitems + i] *= 2.;
    }
  }
  return 0;
}

