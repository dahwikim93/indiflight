#include "cvxgen_ca_cbf_incremental_wrapper.h"
#include "solver.h"

#include <math.h>
#include <stddef.h>

#ifndef ZERO_LIBRARY_MODE
#define CVXGEN_WRAPPER_PRINTF(...) printf(__VA_ARGS__)
#else
#define CVXGEN_WRAPPER_PRINTF(...) do { } while (0)
#endif

static bool initialized = false;

static bool allFinite(const float *values, int count)
{
    for (int i = 0; i < count; i++) {
        if (!isfinite(values[i])) {
            CVXGEN_WRAPPER_PRINTF("cvx err #1");
            return false;
        }
    }

    return true;
}

static double minDouble(double a, double b)
{
    return a < b ? a : b;
}

/*
 * Calculate one original, unscaled CBF value:
 *
 *   -2 rate_cbf' (rateDot_0 - G2*omegaDot_0 + G*du)
 *   + gamma (rateMag_sq - rate_cbf'*rate_cbf)
 *
 * For the split-rate solver, rate_cbf is either rate_tilt=[p,q,0]^T
 * or rate_yaw=[0,0,r]^T.
 */
static double calculateCbfValue(
    const float *rate_cbf,
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
        rateNormSquared += (double)rate_cbf[axis] * rate_cbf[axis];
        rateRateDot += (double)rate_cbf[axis] * rateDot[axis];
    }

    return -2.0 * rateRateDot
        + (double)gamma * ((double)rateMag_sq - rateNormSquared);
}

static double calculateOneCbfMaxMagnitude(
    const float *rate_cbf,
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
        rateNormSquared += (double)rate_cbf[axis] * rate_cbf[axis];
        rateDotDrift += (double)rate_cbf[axis] * drift[axis];
    }

    const double cbfConstant =
        -2.0 * rateDotDrift
        + (double)gamma * ((double)rateMag_sq - rateNormSquared);

    double maximumMagnitude = fabs(cbfConstant);

    for (int actuator = 0; actuator < CVXGEN_CA_CBF_ACTS; actuator++) {
        double coefficient = 0.0;

        for (int axis = 0; axis < CVXGEN_CA_CBF_ATT; axis++) {
            const int index = axis + CVXGEN_CA_CBF_ATT * actuator;
            coefficient += (double)rate_cbf[axis] * G[index];
        }

        coefficient *= -2.0;

        if (fabs(coefficient) > maximumMagnitude) {
            maximumMagnitude = fabs(coefficient);
        }
    }

    return maximumMagnitude;
}

/*
 * Calculate one positive scale for both CBF inequalities.
 *
 * The generated CVXGEN solver shares rateDot_0, G2 and G between the tilt
 * and yaw CBF rows. Therefore, the wrapper scales the two CBF constraints
 * by the same positive factor. This keeps the feasible set unchanged.
 */
static double calculateCbfScale(
    const float *rate_tilt,
    const float *rate_yaw,
    const float *rateDot_0,
    const float *G2,
    const float *omegaDot_0,
    const float *G,
    float gamma_tilt,
    float gamma_yaw,
    float tiltRateMag_sq,
    float yawRateMag_sq)
{
    double maximumMagnitude = calculateOneCbfMaxMagnitude(
        rate_tilt,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        gamma_tilt,
        tiltRateMag_sq
    );

    const double yawMaximumMagnitude = calculateOneCbfMaxMagnitude(
        rate_yaw,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        gamma_yaw,
        yawRateMag_sq
    );

    if (yawMaximumMagnitude > maximumMagnitude) {
        maximumMagnitude = yawMaximumMagnitude;
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
    cvxgenCaCbfInfo_t *info)
{
    if (A == NULL || b == NULL || du_min == NULL || du_max == NULL
        || rate_tilt == NULL || rate_yaw == NULL || rateDot_0 == NULL
        || G2 == NULL || omegaDot_0 == NULL || G == NULL
        || du_out == NULL) {
        CVXGEN_WRAPPER_PRINTF("cvx err #2");
        return false;
    }

    if (!allFinite(A, CVXGEN_CA_CBF_ROWS * CVXGEN_CA_CBF_ACTS)
        || !allFinite(b, CVXGEN_CA_CBF_ROWS)
        || !allFinite(du_min, CVXGEN_CA_CBF_ACTS)
        || !allFinite(du_max, CVXGEN_CA_CBF_ACTS)
        || !allFinite(rate_tilt, CVXGEN_CA_CBF_ATT)
        || !allFinite(rate_yaw, CVXGEN_CA_CBF_ATT)
        || !allFinite(rateDot_0, CVXGEN_CA_CBF_ATT)
        || !allFinite(G2, CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS)
        || !allFinite(omegaDot_0, CVXGEN_CA_CBF_ACTS)
        || !allFinite(G, CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS)
        || !isfinite(gamma_tilt)
        || !isfinite(gamma_yaw)
        || !isfinite(tiltRateMag_sq)
        || !isfinite(yawRateMag_sq)) {
        CVXGEN_WRAPPER_PRINTF("cvx err #3");
        return false;
    }

    if (gamma_tilt < 0.0f || gamma_yaw < 0.0f
        || tiltRateMag_sq < 0.0f || yawRateMag_sq < 0.0f) {
        CVXGEN_WRAPPER_PRINTF("cvx err #4");
        return false;
    }

    for (int i = 0; i < CVXGEN_CA_CBF_ACTS; i++) {
        if (du_min[i] > du_max[i]) {
            CVXGEN_WRAPPER_PRINTF("cvx err #5");
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
        settings.eps = 1e-4;
        settings.resid_tol = 1e-6;
        settings.kkt_reg = 1e-7;

        initialized = true;
    }

    /*
     * Scale A and b together.
     *
     * This changes the objective by one positive constant factor and
     * therefore does not change its minimizer.
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
     * Scale both CBF inequalities by one common positive factor.
     */
    const double cbfScale = calculateCbfScale(
        rate_tilt,
        rate_yaw,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        gamma_tilt,
        gamma_yaw,
        tiltRateMag_sq,
        yawRateMag_sq
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
     * The CBF rate selector vectors remain in their original units.
     *
     * rateDot_0 is scaled as part of the complete CBF inequalities.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ATT; i++) {
        params.rate_tilt[i] = (double)rate_tilt[i];
        params.rate_yaw[i] = (double)rate_yaw[i];
        params.rateDot_0[i] = cbfScale * (double)rateDot_0[i];
    }

    /*
     * G2 and G are scaled by the same common CBF scale.
     */
    for (int i = 0; i < CVXGEN_CA_CBF_ATT * CVXGEN_CA_CBF_ACTS; i++) {
        params.G2[i] = cbfScale * (double)G2[i];
        params.G[i] = cbfScale * (double)G[i];
    }

    /*
     * gamma values are scaled, but rate magnitude limits remain unchanged.
     */
    params.gamma_tilt[0] = cbfScale * (double)gamma_tilt;
    params.gamma_yaw[0] = cbfScale * (double)gamma_yaw;
    params.tiltRateMag_sq[0] = (double)tiltRateMag_sq;
    params.yawRateMag_sq[0] = (double)yawRateMag_sq;

#ifndef ZERO_LIBRARY_MODE
    CVXGEN_WRAPPER_PRINTF(
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
        info->cbf_tilt_value = NAN;
        info->cbf_yaw_value = NAN;
        info->cbf_min_value = NAN;
    }

    if (work.converged != 1) {
        CVXGEN_WRAPPER_PRINTF("cvx err #6");
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
            CVXGEN_WRAPPER_PRINTF("cvx err #7");
            return false;
        }
    }

    /*
     * Calculate the CBF values using the original, unscaled parameters.
     */
    const double tiltCbfValue = calculateCbfValue(
        rate_tilt,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        du_out,
        gamma_tilt,
        tiltRateMag_sq
    );

    const double yawCbfValue = calculateCbfValue(
        rate_yaw,
        rateDot_0,
        G2,
        omegaDot_0,
        G,
        du_out,
        gamma_yaw,
        yawRateMag_sq
    );

    if (info != NULL) {
        info->cbf_tilt_value = tiltCbfValue;
        info->cbf_yaw_value = yawCbfValue;
        info->cbf_min_value = minDouble(tiltCbfValue, yawCbfValue);
    }

    /*
     * Reject a solution that violates either original CBF inequality.
     */
    if (!isfinite(tiltCbfValue) || tiltCbfValue < -1e-5f) {
        CVXGEN_WRAPPER_PRINTF("cvx err #8 tilt, %f", tiltCbfValue);
        return false;
    }

    if (!isfinite(yawCbfValue) || yawCbfValue < -1e-5f) {
        CVXGEN_WRAPPER_PRINTF("cvx err #8 yaw, %f", yawCbfValue);
        return false;
    }

    return true;
}
