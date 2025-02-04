#include "param.h"
#include "./compute_duy.h"
#include "./compute_duy/advx.h"
#include "./compute_duy/advy.h"
#include "./compute_duy/difx.h"
#include "./compute_duy/dify.h"
#include "./compute_duy/pres.h"

int compute_duy(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
#pragma omp parallel for
  for (size_t j = uy_jmin; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      duy[j][i] = 0.;
    }
  }
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  double ** const  p = flow_field-> p;
  const double c = 1. / Re;
  uy_advx(domain, ux, uy, dt, duy);
  uy_advy(domain,     uy, dt, duy);
  uy_difx(domain,  c, uy, dt, duy);
  uy_dify(domain,  c, uy, dt, duy);
  uy_pres(domain,      p, dt, duy);
  return 0;
}

