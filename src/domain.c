#include "domain.h"

const size_t ux_imin = X_PERIODIC ? 1 : 2;
const size_t uy_jmin = Y_PERIODIC ? 1 : 2;

int domain_init(
    domain_t * const domain
) {
  const double lx = 1.;
  const double ly = 3.;
  const size_t nx = 128;
  const size_t ny = 384;
  const double dx = lx / nx;
  const double dy = ly / ny;
  domain->lx = lx;
  domain->ly = ly;
  domain->nx = nx;
  domain->ny = ny;
  domain->dx = dx;
  domain->dy = dy;
  return 0;
}

