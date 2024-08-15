#if !defined(RDFT_H)
#define RDFT_H

#include <stddef.h> // size_t

// planner
typedef struct rdft_plan_t rdft_plan_t;

// create a plan
extern int rdft_init_plan (
    const size_t nitems,
    const size_t repeat_for,
    rdft_plan_t ** const plan
);

// clean-up a plan
extern int rdft_destroy_plan (
    rdft_plan_t ** const plan
);

// perform forward transform
// in : x[0], x[1], ..., x[n-2], x[n-1]
// out: Re(X[0]), Re(X[1]), ..., Re(X[n/2-1]), Re(X[n/2]),
//      Im(X[n/2-1]), Im(X[n/2-2]), ..., Im(X[2]), Im(X[1])
extern int rdft_exec_f (
    rdft_plan_t * const plan,
    double * const xs
);

// perform backward transform
// in : Re(X[0]), Re(X[1]), ..., Re(X[n/2-1]), Re(X[n/2]),
//      Im(X[n/2-1]), Im(X[n/2-2]), ..., Im(X[2]), Im(X[1])
// out: x[0], x[1], ..., x[n-2], x[n-1]
extern int rdft_exec_b (
    rdft_plan_t * const plan,
    double * const xs
);

#endif // RDFT_H
