#if !defined(FLOW_FIELD_H)
#define FLOW_FIELD_H

#include "array.h" // array_t

typedef struct {
  array_t * ux;
  array_t * uy;
  array_t *  p;
  // penalty to enforce zero velocity
  array_t * weight;
} flow_field_t;

extern int flow_field_init (
    flow_field_t * const flow_field
);

extern int flow_field_finalise (
    flow_field_t * const flow_field
);

#endif // FLOW_FIELD_H
