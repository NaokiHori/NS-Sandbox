#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t
#include <stdbool.h> // true, false

#define X_PERIODIC true
#define Y_PERIODIC false

extern const size_t ux_imin;
extern const size_t uy_jmin;

typedef struct {
  double lx;
  double ly;
  size_t nx;
  size_t ny;
  double dx;
  double dy;
} domain_t;

extern int domain_init(
    domain_t * const domain
);

#endif // DOMAIN_H
