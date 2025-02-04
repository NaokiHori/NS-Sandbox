#include <stddef.h> // size_t
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "./integrate.h"
#include "./monitor.h"
#include "./save.h"

typedef struct {
  double monitor;
  double save;
} schedule_t;

int main(
    void
) {
  domain_t domain = {};
  flow_field_t flow_field = {};
  flow_solver_t flow_solver = {};
  if (0 != domain_init(&domain)) {
    return 1;
  }
  if (0 != flow_field_init(&domain, &flow_field)) {
    return 1;
  }
  if (0 != flow_solver_init(&domain, &flow_solver)) {
    return 1;
  }
  const double time_max = 5.e+0;
  const schedule_t rate = {
    .monitor = 1.e-1,
    .save = 2.e-1,
  };
  schedule_t next = rate;
  for (double time = 0.; time < time_max; ) {
    static size_t step = 0;
    double dt = 0.;
    if (0 != integrate(&domain, &flow_field, &flow_solver, &dt)) {
      break;
    }
    step += 1;
    time += dt;
    if (next.monitor < time) {
      monitor(step, time, dt, &domain, &flow_field);
      next.monitor += rate.monitor;
    }
    if (next.save < time) {
      static size_t id = 0;
      save(id, step, time, &domain, &flow_field);
      id += 1;
      next.save += rate.save;
    }
  }
  if (0 != flow_field_finalize(&flow_field)) {
    return 1;
  }
  if (0 != flow_solver_finalize(&flow_solver)) {
    return 1;
  }
  return 0;
}

