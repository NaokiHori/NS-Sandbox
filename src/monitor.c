#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "./monitor.h"

#define ROOT_DIRECTORY "output/log/"

static int output (
    const size_t step,
    const double time,
    const char file_name[],
    const size_t nitems,
    const double * quantities
) {
  errno = 0;
  FILE * const fp = fopen(file_name, "a");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  fprintf(fp, "%10zu % .15e ", step, time);
  for (size_t n = 0; n < nitems; n++) {
    fprintf(fp, "% .15e%c", quantities[n], nitems - 1 == n ? '\n' : ' ');
  }
  fclose(fp);
  return 0;
}

static int print (
    const size_t step,
    const double time,
    const double dt
) {
  FILE * const stream = stdout;
  fprintf(
      stream,
      "step %10zu time % .2e dt % .2e\n",
      step, time, dt
  );
  return 0;
}

static int monitor_divergence (
    const size_t step,
    const double time,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "divergence.dat";
  const array_t * const ux = flow_field->ux;
  const array_t * const uy = flow_field->uy;
  double div_max = 0.;
  double div_sum = 0.;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double dux = - ux[j    ][i    ]
                         + ux[j    ][i + 1];
      const double duy = - uy[j    ][i    ]
                         + uy[j + 1][i    ];
      const double div = (
          + 1. / DX * dux
          + 1. / DY * duy
      );
      div_max = fmax(div_max, fabs(div));
      div_sum = div_sum + div;
    }
  }
  return output(step, time, file_name, 2, (double []){div_max, div_sum});
}

static int monitor_max_velocity (
    const size_t step,
    const double time,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "max_velocity.dat";
  const array_t * const ux = flow_field->ux;
  const array_t * const uy = flow_field->uy;
  double ux_max = 0.;
  double uy_max = 0.;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      ux_max = fmax(ux_max, fabs(ux[j    ][i    ]));
    }
  }
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      uy_max = fmax(uy_max, fabs(uy[j    ][i    ]));
    }
  }
  return output(step, time, file_name, 2, (double []){ux_max, uy_max});
}

int monitor (
    const size_t step,
    const double time,
    const double dt,
    const flow_field_t * const flow_field
) {
  if (0 != print(step, time, dt)) {
    LOGGER_FAILURE("failed to output metrics");
    goto abort;
  }
  if (0 != monitor_divergence(step, time, flow_field)) {
    LOGGER_FAILURE("failed to check / output divergence");
    goto abort;
  }
  if (0 != monitor_max_velocity(step, time, flow_field)) {
    LOGGER_FAILURE("failed to check / output maximum velocity");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to monitor flow field");
  return 1;
}

