#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "dft/rdft.h"

static const double pi = 3.14159265358979324;

struct rdft_plan_t {
  // length of real signal
  // NOTE: for inverse transform, the input complex signal
  //       should have nitems / 2 + 1 elements
  size_t nitems;
  // repeat DFTs for specified times
  size_t repeat_for;
  // pre-computed cosine / sine values
  double * table_cos;
  double * table_sin;
  // internal buffer
  double complex * buf;
};

static void * memory_alloc(
    const size_t size
) {
  void * const ptr = malloc(size);
  if (NULL == ptr) {
    fprintf(stderr, "[FATAL %s:%d] failed to allocate %zu bytes\n", __FILE__, __LINE__, size);
    return NULL;
  }
  return ptr;
}

static void memory_free(
    void * const ptr
) {
  free(ptr);
}

// recursive Cooley-Tukey FFT for complex input/output
static int dft(
    const size_t nitems,
    const double sign,
    const size_t stride,
    const double * const table_cos,
    const double * const table_sin,
    const double complex * xs,
    double complex * ys
) {
  if (1 == nitems) {
    ys[0] = xs[0];
  } else if (0 == nitems % 2) {
    dft(nitems / 2, sign, stride * 2, table_cos, table_sin, xs         , ys             );
    dft(nitems / 2, sign, stride * 2, table_cos, table_sin, xs + stride, ys + nitems / 2);
    for (size_t i = 0; i < nitems / 2; i++) {
      const size_t j = i + nitems / 2;
      const double c = table_cos[2 * stride * i];
      const double s = table_sin[2 * stride * i];
      const double complex twiddle = c + sign * I * s;
      const double complex e = ys[i];
      const double complex o = ys[j] * twiddle;
      ys[i] = e + o;
      ys[j] = e - o;
    }
  } else {
    // naive O(N^2) DFT
    for (size_t k = 0; k < nitems; k++) {
      double complex * y = ys + k;
      *y = 0. + I * 0.;
      for (size_t n = 0; n < nitems; n++) {
        *y += xs[stride * n] * cexp(sign * 2. * pi * n * k * I / nitems);
      }
    }
  }
	return 0;
}

int rdft_exec_f(
    rdft_plan_t * const plan,
    double * const xs
) {
  if (NULL == plan) {
    puts("uninitialized plan is passed");
    return 1;
  }
  const size_t nitems = plan->nitems;
  const size_t repeat_for = plan->repeat_for;
  const double * const table_cos = plan->table_cos;
  const double * const table_sin = plan->table_sin;
  double complex * const zs = plan->buf;
#pragma omp parallel for
  for (size_t j = 0; j < repeat_for; j++) {
    double * xs_j = xs + j * nitems;
    double complex * zs_j = zs + j * (nitems / 2 + 1);
    // create a signal composed of N/2 complex numbers
    // x[2n] + I x[2n + 1] (n = 0, 1, ..., N / 2 - 1)
    // NOTE: the original memory layout already satisfies the requirement
    //       due to C99 standard, so we just cast and use it
    dft(nitems / 2, - 1., 1, table_cos, table_sin, (double complex *)xs_j, zs_j);
    // duplicate for later convenience
    zs_j[nitems / 2] = zs_j[0];
    // from the fourier transformed signal, compute FFT of even / odd signals
    for (size_t i = 0; i < nitems / 2 + 1; i++) {
      const double complex e = + 0.5 * zs_j[i] + 0.5 * conj(zs_j[nitems / 2 - i]);
      const double complex o = - 0.5 * zs_j[i] + 0.5 * conj(zs_j[nitems / 2 - i]);
      const double c = table_cos[i];
      const double s = table_sin[i];
      const double complex twiddle = c - I * s;
      const double complex result = e + o * I * twiddle;
      xs_j[i] = creal(result);
      if (0 != i && nitems / 2 != i) {
        xs_j[nitems - i] = cimag(result);
      }
    }
  }
  return 0;
}

int rdft_exec_b(
    rdft_plan_t * const plan,
    double * const xs
) {
  if (NULL == plan) {
    puts("uninitialized plan is passed");
    return 1;
  }
  const size_t nitems = plan->nitems;
  const size_t repeat_for = plan->repeat_for;
  const double * const table_cos = plan->table_cos;
  const double * const table_sin = plan->table_sin;
  double complex * const zs = plan->buf;
#pragma omp parallel for
  for (size_t j = 0; j < repeat_for; j++) {
    double * xs_j = xs + j * nitems;
    double complex * zs_j = zs + j * (nitems / 2 + 1);
    for (size_t i = 0; i < nitems / 2; i++) {
      const double real0 =               xs_j[             i];
      const double imag0 = 0 == i ? 0. : xs_j[nitems     - i];
      const double real1 =               xs_j[nitems / 2 - i];
      const double imag1 = 0 == i ? 0. : xs_j[nitems / 2 + i];
      const double complex val0 = real0 + I * imag0;
      const double complex val1 = real1 + I * imag1;
      const double complex e = + 0.5 * val0 + 0.5 * conj(val1);
      const double complex o = + 0.5 * val0 - 0.5 * conj(val1);
      const double c = table_cos[i];
      const double s = table_sin[i];
      const double complex twiddle = c + I * s;
      zs_j[i] = e + o * I * twiddle;
    }
    dft(nitems / 2, + 1., 1, table_cos, table_sin, zs_j, (double complex *)xs_j);
    // NOTE: performing DFTs whose size is nitems / 2
    //       halves the amplitude of the resulting signal,
    //       which is compensated here
    for (size_t i = 0; i < nitems; i++) {
      xs_j[i] *= 2.;
    }
  }
  return 0;
}

int rdft_init_plan(
    const size_t nitems,
    const size_t repeat_for,
    rdft_plan_t ** const plan
) {
  if (0 != nitems % 2) {
    printf("signal length (%zu) should be a multiple of 2\n", nitems);
    return 1;
  }
  *plan = memory_alloc(1 * sizeof(rdft_plan_t));
  (*plan)->nitems = nitems;
  (*plan)->repeat_for = repeat_for;
  double ** table_cos = &(*plan)->table_cos;
  double ** table_sin = &(*plan)->table_sin;
  double complex ** buf = &(*plan)->buf;
  *table_cos = memory_alloc((nitems / 2 + 1) * sizeof(double));
  *table_sin = memory_alloc((nitems / 2 + 1) * sizeof(double));
  *buf = memory_alloc((nitems / 2 + 1) * repeat_for * sizeof(double complex));
  // prepare cosine / sine tables
  for (size_t i = 0; i < nitems / 2 + 1; i++) {
    (*table_cos)[i] = cos(2. * pi * i / nitems);
    (*table_sin)[i] = sin(2. * pi * i / nitems);
  }
  return 0;
}

// clean-up a plan
int rdft_destroy_plan(
    rdft_plan_t ** const plan
) {
  memory_free((*plan)->table_cos);
  memory_free((*plan)->table_sin);
  memory_free((*plan)->buf);
  memory_free(*plan);
  *plan = NULL;
  return 0;
}

