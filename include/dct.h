#if !defined(DCT_H)
#define DCT_H

#include <stddef.h> // size_t

// planner
typedef struct dct_plan_t dct_plan_t;

// create a plan
extern int dct_init_plan (
    const size_t nitems,
    const size_t repeat_for,
    dct_plan_t ** plan
);

// clean-up a plan
extern int dct_destroy_plan (
    dct_plan_t ** plan
);

// perform forward transform (DCT type 2)
extern int dct_exec_f (
    dct_plan_t * plan,
    double * restrict xs
);

// perform backward transform (DCT type 3)
extern int dct_exec_b (
    dct_plan_t * plan,
    double * restrict xs
);

#endif // DCT_H
