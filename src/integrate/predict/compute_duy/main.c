#include "param.h"
#include "../compute_duy.h"
#include "./advx.h"
#include "./advy.h"
#include "./difx.h"
#include "./dify.h"
#include "./pres.h"

int compute_duy (
    const flow_field_t * const flow_field,
    const double dt,
    array_t * const duy
) {
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      duy[j][i] = 0.;
    }
  }
  const array_t * const ux = flow_field->ux;
  const array_t * const uy = flow_field->uy;
  const array_t * const  p = flow_field-> p;
  const double c = 1. / Re;
  uy_advx(ux, uy, dt, duy);
  uy_advy(    uy, dt, duy);
  uy_difx( c, uy, dt, duy);
  uy_dify( c, uy, dt, duy);
  uy_pres(     p, dt, duy);
  return 0;
}

