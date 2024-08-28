#include "param.h"
#include "../compute_dux.h"
#include "./advx.h"
#include "./advy.h"
#include "./difx.h"
#include "./dify.h"
#include "./pres.h"

int compute_dux (
    const flow_field_t * const flow_field,
    const double dt,
    array_t * const dux
) {
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = ux_imin; i <= NX; i++) {
      dux[j][i] = 0.;
    }
  }
  const array_t * const ux = flow_field->ux;
  const array_t * const uy = flow_field->uy;
  const array_t * const  p = flow_field-> p;
  const double c = 1. / Re;
  ux_advx(    ux, dt, dux);
  ux_advy(uy, ux, dt, dux);
  ux_difx( c, ux, dt, dux);
  ux_dify( c, ux, dt, dux);
  ux_pres(     p, dt, dux);
  return 0;
}

