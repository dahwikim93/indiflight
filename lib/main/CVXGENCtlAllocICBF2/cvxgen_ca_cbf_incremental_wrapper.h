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
    double cbf_tilt_value;
    double cbf_yaw_value;
    double cbf_min_value;
} cvxgenCaCbfInfo_t;

/*
 * Solves:
 *
 *   minimize    ||A*du - b||^2
 *
 *   subject to  du_min <= du <= du_max
 *
 *               -2*rate_tilt'*(rateDot_0 - G2*omegaDot_0 + G*du)
 *               + gamma_tilt*(tiltRateMag_sq - rate_tilt'*rate_tilt) >= 0
 *
 *               -2*rate_yaw'*(rateDot_0 - G2*omegaDot_0 + G*du)
 *               + gamma_yaw*(yawRateMag_sq - rate_yaw'*rate_yaw) >= 0
 *
 * Matrix storage:
 *   A, G and G2 must be column-major.
 *
 * Dimensions:
 *   A                  10x4
 *   b                  10x1
 *   du_min              4x1
 *   du_max              4x1
 *   rate_tilt           3x1   usually [p, q, 0]^T
 *   rate_yaw            3x1   usually [0, 0, r]^T
 *   rateDot_0           3x1
 *   G2                  3x4
 *   omegaDot_0          4x1
 *   G                   3x4
 *   gamma_tilt          scalar
 *   gamma_yaw           scalar
 *   tiltRateMag_sq      scalar, e.g. tilt_rate_max^2 in rad^2/s^2
 *   yawRateMag_sq       scalar, e.g. yaw_rate_max^2 in rad^2/s^2
 */
bool cvxgenCaCbfSolve(
    const float *A,
    const float *b,
    const float *du_min,
    const float *du_max,
    const float *rate_tilt,
    const float *rate_yaw,
    const float *rateDot_0,
    const float *G2,
    const float *omegaDot_0,
    const float *G,
    float gamma_tilt,
    float gamma_yaw,
    float tiltRateMag_sq,
    float yawRateMag_sq,
    float *du_out,
    cvxgenCaCbfInfo_t *info
);

#endif
