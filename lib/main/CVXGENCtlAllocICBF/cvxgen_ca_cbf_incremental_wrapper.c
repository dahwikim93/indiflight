#include "cvxgen_ca_cbf_incremental_wrapper.h"
#include "solver.h"

#include <math.h>
#include <stddef.h>

static bool initialized = false;

static bool allFinite(const float *values, int count)
{
    for (int i = 0; i < count; i++) {
        if (!isfinite(values[i])) {
            printf("cvx err #1");
            return false;
        }
    }

    return true;
}

/*
 * Calculate the original, unscaled CBF value:
 *
 *   -2 rate' (rateDot_0 - G2*omegaDot_0 + G*du)
 *   + gamma (rateMag_sq - rate'*rate)
 */
static double calculateCbfValue(
    const float *rate,
    const float *rateDot_0,
    const float *G2,
    const float *omegaDot_0,
    const float *G,
    const float *du,
    float gamma,
    float rateMag_sq)
{
    double rateDot[CVXGEN_CA_CBF_ATT];

    for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
        double value = (double)rateDot_0[axis];

        for (int actuator = 0; actuator < CVXGEN_CA_CBF_ACTS; actuator++) {
            const int index = axis + CVXGEN_CA_CBF_ATT * actuator;

            value -= (double)G2[index] * omegaDot_0[actuator];
            value += (double)G[index] * du[actuator];
        }

        rateDot[axis] = value;
    }

    double rateNormSquared = 0.0;
    double rateRateDot = 0.0;

    for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
        rateNormSquared += (double)rate[axis] * rate[axis];
        rateRateDot += (double)rate[axis] * rateDot[axis];
    }

    return -2.0 * rateRateDot
        + (double)gamma * ((double)rateMag_sq - rateNormSquared);
}

/*
 * Calculate a positive scale for the complete CBF inequality.
 *
 * Original CBF:
 *
 *   c + a' du >= 0
 *
 * where:
 *
 *   c = -2 rate' (rateDot_0 - G2*omegaDot_0)
 *       + gamma (rateMag_sq - rate'*rate)
 *
 *   a_i = -2 rate' G(:,i)
 *
 * rateDot_0, G2, G and gamma are all multiplied by this same
 * positive scale. Therefore, the feasible set is unchanged.
 */
static double calculateCbfScale(
    const float *rate,
    const float *rateDot_0,
    const float *G2,
    const float *omegaDot_0,
    const float *G,
    float gamma,
    float rateMag_sq)
{
    double drift[CVXGEN_CA_CBF_ATT];

    for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
        drift[axis] = (double)rateDot_0[axis];

        for (int actuator = 0; actuator < CVXGEN_CA_CBF_ACTS; actuator++) {
            const int index = axis + CVXGEN_CA_CBF_ATT * actuator;
            drift[axis] -= (double)G2[index] * omegaDot_0[actuator];
        }
    }

    double rateNormSquared = 0.0;
    double rateDotDrift = 0.0;

    for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
        rateNormSquared += (double)rate[axis] * rate[axis];
        rateDotDrift += (double)rate[axis] * drift[axis];
    }

    const double cbfConstant =
        -2.0 * rateDotDrift
        + (double)gamma * ((double)rateMag_sq - rateNormSquared);

    double maximumMagnitude = fabs(cbfConstant);

    for (int actuator = 0; actuator < CVXGEN_CA_CBF_ACTS; actuator++) {
        double coefficient = 0.0;

        for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
            const int index = axis + CVXGEN_CA_CBF_ATT * actuator;
            coefficient += (double)rate[axis] * G[index];
        }

        coefficient *= -2.0;

        if (fabs(coefficient) > maximumMagnitude) {
            maximumMagnitude = fabs(coefficient);
        }
    }

    /*
     * Do not scale small CBF constraints upward.
     */
    if (maximumMagnitude <= 1.0) {
        return 1.0;
    }

    return 1.0 / maximumMagnitude;
}

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
    cvxgenCaCbfInfo_t *info)
{
    if (A == NULL || b == NULL || du_min == NULL || du_max == NULL
        || rate == NULL || rateDot_0 == NULL || G2 == NULL
        || omegaDot_0 == NULL || G == NULL || du_out == NULL) {
        printf("cvx err #2");
        return false;
    }

    if (!allFinite(A, CVXGEN_CA_CBF_ROWS * CVXGEN_CA_CBF_ACTS)
        || !allFinite(b, CVXGEN_CA_CBF_ROWS)
        || !allFinite(du_min, CVXGEN_CA_CBF_ACTS)
        || !allFinite(du_max, CVXGEN_CA_CBF_ACTS)
        || !allFinite(rate, CVXGEN_CA_CBF_ATT)
        || !allFinite(rateDot_0, CVXGEN_CA_CBF_ATT)
        || !allFinite(G2, CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS)
        || !allFinite(omegaDot_0, CVXGEN_CA_CBF_ACTS)
        || !allFinite(G, CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS)
        || !isfinite(gamma)
        || !isfinite(rateMag_sq)) {
        printf("cvx err #3");
        return false;
    }

    if (gamma < 0.0f || rateMag_sq < 0.0f) {
        printf("cvx err #4");
        return false;
    }

    for (int i = 0; i < CVXGEN_CA_CBF_ACTS; i++) {
        if (du_min[i] > du_max[i]) {
            printf("cvx err #5");
            return false;
        }
    }

    if (!initialized) {
        set_defaults();
        setup_indexing();

        settings.verbose = 1;
        settings.max_iters = 30;
        settings.refine_steps = 1;
        settings.better_start = 1;
        settings.eps = 1e-6;
        settings.resid_tol = 1e-6;
        settings.kkt_reg = 1e-7;

        initialized = true;
    }

    /*
     * Scale A and b together.
     *
     * This changes the objective by one positive constant factor
     * and therefore does not change its minimizer.
     */
    double objectiveScale = 1.0;

    for (int i = 0; i < CVXGEN_CA_CBF_ROWS * CVXGEN_CA_CBF_ACTS; i++) {
        const double magnitude = fabs((double)A[i]);

        if (magnitude > objectiveScale) {
            objectiveScale = magnitude;
        }
    }

    for (int i = 0; i < CVXGEN_CA_CBF_ROWS; i++) {
        const double magnitude = fabs((double)b[i]);

        if (magnitude > objectiveScale) {
            objectiveScale = magnitude;
        }
    }

    /*
     * Scale the complete CBF inequality by one positive factor.
     */
    const double cbfScale = calculateCbfScale(
        rate,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        gamma,
        rateMag_sq
    );

    /*
     * Objective parameters.
     *
     * CVXGEN dense matrices use column-major storage.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ROWS * CVXGEN_CA_CBF_ACTS; i++) {
        params.A[i] = (double)A[i] / objectiveScale;
    }

    for (int i = 0; i < CVXGEN_CA_CBF_ROWS; i++) {
        params.b[i] = (double)b[i] / objectiveScale;
    }

    /*
     * Actuator bounds and current rotor acceleration are not scaled.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ACTS; i++) {
        params.du_min[i] = (double)du_min[i];
        params.du_max[i] = (double)du_max[i];
        params.omegaDot_0[i] = (double)omegaDot_0[i];
    }

    /*
     * rate remains in its original units.
     *
     * rateDot_0 is scaled as part of the complete CBF inequality.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ATT; i++) {
        params.rate[i] = (double)rate[i];
        params.rateDot_0[i] = cbfScale * (double)rateDot_0[i];
    }

    /*
     * G2 and G are scaled by the same CBF scale.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS; i++) {
        params.G2[i] = cbfScale * (double)G2[i];
        params.G[i] = cbfScale * (double)G[i];
    }

    /*
     * gamma is scaled, but rateMag_sq remains unchanged.
     */
    params.gamma[0] = cbfScale * (double)gamma;
    params.rateMag_sq[0] = (double)rateMag_sq;

#ifndef ZERO_LIBRARY_MODE
    printf(
        "CVXGEN scales: objective=% .6e cbf=% .6e\n",
        objectiveScale,
        cbfScale
    );
#endif

    const long iterations = solve();

    if (info != NULL) {
        info->iterations = (int)iterations;
        info->gap = work.gap;
        info->inequality_residual_squared = work.ineq_resid_squared;
        info->objective_value = work.optval;
        info->cbf_value = NAN;
    }

    if (work.converged != 1) {
        printf("cvx err #6");
        return false;
    }

    for (int i = 0; i < CVXGEN_CA_CBF_ACTS; i++) {
        du_out[i] = (float)vars.du_cvx[i];
    }

    /*
     * Check bounds explicitly after convergence.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ACTS; i++) {
        if (du_out[i] < du_min[i] - 1e-5f
            || du_out[i] > du_max[i] + 1e-5f) {
            printf("cvx err #7");
            return false;
        }
    }

    /*
     * Calculate the CBF value using the original, unscaled parameters.
     */
    const double cbfValue = calculateCbfValue(
        rate,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        du_out,
        gamma,
        rateMag_sq
    );

    if (info != NULL) {
        info->cbf_value = cbfValue;
    }

    /*
     * Reject a solution that violates the original CBF inequality.
     */
    if (!isfinite(cbfValue) || cbfValue < -1e-3) {
        printf("cvx err #8, %f", cbfValue);
        return false;
    }

    return true;
}
