#pragma once

#include <stdbool.h>

#define CVXGEN_CA_ROWS 10
#define CVXGEN_CA_ACTS 4

bool cvxgenControlAllocationSolve(
    const float *A_as,
    const float *b_as,
    const float *du_min,
    const float *du_max,
    float *du_out,
    int *num_iters,
    double *gap,
    double *ineq_resid_squared
);