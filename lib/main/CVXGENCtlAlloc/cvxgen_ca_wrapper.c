#include "cvxgen_ca_wrapper.h"
#include "solver.h"

bool cvxgenControlAllocationSolve(
    const float *A_as,
    const float *b_as,
    const float *du_min,
    const float *du_max,
    float *du_out,
    int *num_iters,
    double *gap,
    double *ineq_resid_squared
) {
    static bool initialized = false;

    if (!initialized) {
        set_defaults();
        setup_indexing();

        settings.verbose = 0;
        settings.max_iters = 30;
        settings.eps = 1e-4;
        settings.refine_steps = 0;

        initialized = true;
    }

    for (int i = 0; i < CVXGEN_CA_ROWS * CVXGEN_CA_ACTS; i++) {
        params.A[i] = (double)A_as[i];
    }

    for (int i = 0; i < CVXGEN_CA_ROWS; i++) {
        params.b[i] = (double)b_as[i];
    }

    for (int i = 0; i < CVXGEN_CA_ACTS; i++) {
        params.du_min[i] = (double)du_min[i];
        params.du_max[i] = (double)du_max[i];
    }

    int iters = solve();

    if (num_iters != 0) {
        *num_iters = iters;
    }

    if (gap != 0) {
        *gap = work.gap;
    }

    if (ineq_resid_squared != 0) {
        *ineq_resid_squared = work.ineq_resid_squared;
    }

    if (work.converged != 1) {
        return false;
    }

    for (int i = 0; i < CVXGEN_CA_ACTS; i++) {
        du_out[i] = (float)vars.du_cvx[i];
    }

    return true;
}