#ifndef CVXGEN_CA_CBF_INCREMENTAL_WRAPPER_H
#define CVXGEN_CA_CBF_INCREMENTAL_WRAPPER_H

#include <stdbool.h>

#define CVXGEN_CA_CBF_ROWS 10
#define CVXGEN_CA_CBF_ACTS 4
#define CVXGEN_CA_CBF_ATT 3

typedef struct {
    int iterations;
    double gap;
    double inequality_residual_squared;
    double objective_value;
    double cbf_value;
} cvxgenCaCbfInfo_t;

/*
 * Solves:
 *
 *   minimize    ||A*du - b||^2
 *
 *   subject to  du_min <= du <= du_max
 *
 *               -2*rate'*(rateDot_0 - G2*omegaDot_0 + G*du)
 *               + gamma*(rateMag_sq - rate'*rate) >= 0
 *
 * Matrix storage:
 *   A, G and G2 must be column-major.
 *
 * Dimensions:
 *   A           10x4
 *   b           10x1
 *   du_min       4x1
 *   du_max       4x1
 *   rate         3x1
 *   rateDot_0    3x1
 *   G2           3x4
 *   omegaDot_0   4x1
 *   G            3x4
 */
bool cvxgenCaCbfSolve(
    const float *A,
    const float *b,
    const float *du_min,
    const float *du_max,
    const float *rate,
    const float *rateDot_0,
    const float *G2,
    const float *omegaDot_0,
    const float *G,
    float gamma,
    float rateMag_sq,
    float *du_out,
    cvxgenCaCbfInfo_t *info
);

#endif
