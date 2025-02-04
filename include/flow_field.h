#if !defined(FLOW_FIELD_H)
#define FLOW_FIELD_H

#include "domain.h" // domain_t

typedef struct {
  double ** ux;
  double ** uy;
  double **  p;
  // penalty to enforce zero velocity
  double ** weight;
} flow_field_t;

extern int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
);

extern int flow_field_finalize(
    flow_field_t * const flow_field
);

#endif // FLOW_FIELD_H
