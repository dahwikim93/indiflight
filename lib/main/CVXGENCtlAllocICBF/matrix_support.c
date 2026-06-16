/* Produced by CVXGEN, 2026-06-16 09:21:22 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
  (void)lhs;
  (void)rhs;
}
void multbymAT(double *lhs, double *rhs) {
  (void)rhs;

  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[1]*(-1);
  lhs[2] = -rhs[2]*(-1);
  lhs[3] = -rhs[3]*(-1);
  lhs[4] = -rhs[0]*(1);
  lhs[5] = -rhs[1]*(1);
  lhs[6] = -rhs[2]*(1);
  lhs[7] = -rhs[3]*(1);
  lhs[8] = -rhs[0]*(2*(params.rate[0]*params.G[0]+params.rate[1]*params.G[1]+params.rate[2]*params.G[2]))-rhs[1]*(2*(params.rate[0]*params.G[3]+params.rate[1]*params.G[4]+params.rate[2]*params.G[5]))-rhs[2]*(2*(params.rate[0]*params.G[6]+params.rate[1]*params.G[7]+params.rate[2]*params.G[8]))-rhs[3]*(2*(params.rate[0]*params.G[9]+params.rate[1]*params.G[10]+params.rate[2]*params.G[11]));
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[4]*(1)-rhs[8]*(2*(params.rate[0]*params.G[0]+params.rate[1]*params.G[1]+params.rate[2]*params.G[2]));
  lhs[1] = -rhs[1]*(-1)-rhs[5]*(1)-rhs[8]*(2*(params.rate[0]*params.G[3]+params.rate[1]*params.G[4]+params.rate[2]*params.G[5]));
  lhs[2] = -rhs[2]*(-1)-rhs[6]*(1)-rhs[8]*(2*(params.rate[0]*params.G[6]+params.rate[1]*params.G[7]+params.rate[2]*params.G[8]));
  lhs[3] = -rhs[3]*(-1)-rhs[7]*(1)-rhs[8]*(2*(params.rate[0]*params.G[9]+params.rate[1]*params.G[10]+params.rate[2]*params.G[11]));
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*work.quad_601932361728[0])+rhs[1]*(2*work.quad_601932361728[4])+rhs[2]*(2*work.quad_601932361728[8])+rhs[3]*(2*work.quad_601932361728[12]);
  lhs[1] = rhs[0]*(2*work.quad_601932361728[1])+rhs[1]*(2*work.quad_601932361728[5])+rhs[2]*(2*work.quad_601932361728[9])+rhs[3]*(2*work.quad_601932361728[13]);
  lhs[2] = rhs[0]*(2*work.quad_601932361728[2])+rhs[1]*(2*work.quad_601932361728[6])+rhs[2]*(2*work.quad_601932361728[10])+rhs[3]*(2*work.quad_601932361728[14]);
  lhs[3] = rhs[0]*(2*work.quad_601932361728[3])+rhs[1]*(2*work.quad_601932361728[7])+rhs[2]*(2*work.quad_601932361728[11])+rhs[3]*(2*work.quad_601932361728[15]);
}
void fillq(void) {
  work.q[0] = -2*(params.A[0]*params.b[0]+params.A[1]*params.b[1]+params.A[2]*params.b[2]+params.A[3]*params.b[3]+params.A[4]*params.b[4]+params.A[5]*params.b[5]+params.A[6]*params.b[6]+params.A[7]*params.b[7]+params.A[8]*params.b[8]+params.A[9]*params.b[9]);
  work.q[1] = -2*(params.A[10]*params.b[0]+params.A[11]*params.b[1]+params.A[12]*params.b[2]+params.A[13]*params.b[3]+params.A[14]*params.b[4]+params.A[15]*params.b[5]+params.A[16]*params.b[6]+params.A[17]*params.b[7]+params.A[18]*params.b[8]+params.A[19]*params.b[9]);
  work.q[2] = -2*(params.A[20]*params.b[0]+params.A[21]*params.b[1]+params.A[22]*params.b[2]+params.A[23]*params.b[3]+params.A[24]*params.b[4]+params.A[25]*params.b[5]+params.A[26]*params.b[6]+params.A[27]*params.b[7]+params.A[28]*params.b[8]+params.A[29]*params.b[9]);
  work.q[3] = -2*(params.A[30]*params.b[0]+params.A[31]*params.b[1]+params.A[32]*params.b[2]+params.A[33]*params.b[3]+params.A[34]*params.b[4]+params.A[35]*params.b[5]+params.A[36]*params.b[6]+params.A[37]*params.b[7]+params.A[38]*params.b[8]+params.A[39]*params.b[9]);
}
void fillh(void) {
  work.h[0] = -params.du_min[0];
  work.h[1] = -params.du_min[1];
  work.h[2] = -params.du_min[2];
  work.h[3] = -params.du_min[3];
  work.h[4] = params.du_max[0];
  work.h[5] = params.du_max[1];
  work.h[6] = params.du_max[2];
  work.h[7] = params.du_max[3];
  work.h[8] = -2*(params.rate[0]*(params.rateDot_0[0]-(params.G2[0]*params.omegaDot_0[0]+params.G2[3]*params.omegaDot_0[1]+params.G2[6]*params.omegaDot_0[2]+params.G2[9]*params.omegaDot_0[3]))+params.rate[1]*(params.rateDot_0[1]-(params.G2[1]*params.omegaDot_0[0]+params.G2[4]*params.omegaDot_0[1]+params.G2[7]*params.omegaDot_0[2]+params.G2[10]*params.omegaDot_0[3]))+params.rate[2]*(params.rateDot_0[2]-(params.G2[2]*params.omegaDot_0[0]+params.G2[5]*params.omegaDot_0[1]+params.G2[8]*params.omegaDot_0[2]+params.G2[11]*params.omegaDot_0[3])))+params.gamma[0]*(params.rateMag_sq[0]-(params.rate[0]*params.rate[0]+params.rate[1]*params.rate[1]+params.rate[2]*params.rate[2]));
}
void fillb(void) {
}
void pre_ops(void) {
  work.quad_601932361728[0] = params.A[0]*params.A[0]+params.A[1]*params.A[1]+params.A[2]*params.A[2]+params.A[3]*params.A[3]+params.A[4]*params.A[4]+params.A[5]*params.A[5]+params.A[6]*params.A[6]+params.A[7]*params.A[7]+params.A[8]*params.A[8]+params.A[9]*params.A[9];
  work.quad_601932361728[4] = params.A[0]*params.A[10]+params.A[1]*params.A[11]+params.A[2]*params.A[12]+params.A[3]*params.A[13]+params.A[4]*params.A[14]+params.A[5]*params.A[15]+params.A[6]*params.A[16]+params.A[7]*params.A[17]+params.A[8]*params.A[18]+params.A[9]*params.A[19];
  work.quad_601932361728[8] = params.A[0]*params.A[20]+params.A[1]*params.A[21]+params.A[2]*params.A[22]+params.A[3]*params.A[23]+params.A[4]*params.A[24]+params.A[5]*params.A[25]+params.A[6]*params.A[26]+params.A[7]*params.A[27]+params.A[8]*params.A[28]+params.A[9]*params.A[29];
  work.quad_601932361728[12] = params.A[0]*params.A[30]+params.A[1]*params.A[31]+params.A[2]*params.A[32]+params.A[3]*params.A[33]+params.A[4]*params.A[34]+params.A[5]*params.A[35]+params.A[6]*params.A[36]+params.A[7]*params.A[37]+params.A[8]*params.A[38]+params.A[9]*params.A[39];
  work.quad_601932361728[1] = params.A[10]*params.A[0]+params.A[11]*params.A[1]+params.A[12]*params.A[2]+params.A[13]*params.A[3]+params.A[14]*params.A[4]+params.A[15]*params.A[5]+params.A[16]*params.A[6]+params.A[17]*params.A[7]+params.A[18]*params.A[8]+params.A[19]*params.A[9];
  work.quad_601932361728[5] = params.A[10]*params.A[10]+params.A[11]*params.A[11]+params.A[12]*params.A[12]+params.A[13]*params.A[13]+params.A[14]*params.A[14]+params.A[15]*params.A[15]+params.A[16]*params.A[16]+params.A[17]*params.A[17]+params.A[18]*params.A[18]+params.A[19]*params.A[19];
  work.quad_601932361728[9] = params.A[10]*params.A[20]+params.A[11]*params.A[21]+params.A[12]*params.A[22]+params.A[13]*params.A[23]+params.A[14]*params.A[24]+params.A[15]*params.A[25]+params.A[16]*params.A[26]+params.A[17]*params.A[27]+params.A[18]*params.A[28]+params.A[19]*params.A[29];
  work.quad_601932361728[13] = params.A[10]*params.A[30]+params.A[11]*params.A[31]+params.A[12]*params.A[32]+params.A[13]*params.A[33]+params.A[14]*params.A[34]+params.A[15]*params.A[35]+params.A[16]*params.A[36]+params.A[17]*params.A[37]+params.A[18]*params.A[38]+params.A[19]*params.A[39];
  work.quad_601932361728[2] = params.A[20]*params.A[0]+params.A[21]*params.A[1]+params.A[22]*params.A[2]+params.A[23]*params.A[3]+params.A[24]*params.A[4]+params.A[25]*params.A[5]+params.A[26]*params.A[6]+params.A[27]*params.A[7]+params.A[28]*params.A[8]+params.A[29]*params.A[9];
  work.quad_601932361728[6] = params.A[20]*params.A[10]+params.A[21]*params.A[11]+params.A[22]*params.A[12]+params.A[23]*params.A[13]+params.A[24]*params.A[14]+params.A[25]*params.A[15]+params.A[26]*params.A[16]+params.A[27]*params.A[17]+params.A[28]*params.A[18]+params.A[29]*params.A[19];
  work.quad_601932361728[10] = params.A[20]*params.A[20]+params.A[21]*params.A[21]+params.A[22]*params.A[22]+params.A[23]*params.A[23]+params.A[24]*params.A[24]+params.A[25]*params.A[25]+params.A[26]*params.A[26]+params.A[27]*params.A[27]+params.A[28]*params.A[28]+params.A[29]*params.A[29];
  work.quad_601932361728[14] = params.A[20]*params.A[30]+params.A[21]*params.A[31]+params.A[22]*params.A[32]+params.A[23]*params.A[33]+params.A[24]*params.A[34]+params.A[25]*params.A[35]+params.A[26]*params.A[36]+params.A[27]*params.A[37]+params.A[28]*params.A[38]+params.A[29]*params.A[39];
  work.quad_601932361728[3] = params.A[30]*params.A[0]+params.A[31]*params.A[1]+params.A[32]*params.A[2]+params.A[33]*params.A[3]+params.A[34]*params.A[4]+params.A[35]*params.A[5]+params.A[36]*params.A[6]+params.A[37]*params.A[7]+params.A[38]*params.A[8]+params.A[39]*params.A[9];
  work.quad_601932361728[7] = params.A[30]*params.A[10]+params.A[31]*params.A[11]+params.A[32]*params.A[12]+params.A[33]*params.A[13]+params.A[34]*params.A[14]+params.A[35]*params.A[15]+params.A[36]*params.A[16]+params.A[37]*params.A[17]+params.A[38]*params.A[18]+params.A[39]*params.A[19];
  work.quad_601932361728[11] = params.A[30]*params.A[20]+params.A[31]*params.A[21]+params.A[32]*params.A[22]+params.A[33]*params.A[23]+params.A[34]*params.A[24]+params.A[35]*params.A[25]+params.A[36]*params.A[26]+params.A[37]*params.A[27]+params.A[38]*params.A[28]+params.A[39]*params.A[29];
  work.quad_601932361728[15] = params.A[30]*params.A[30]+params.A[31]*params.A[31]+params.A[32]*params.A[32]+params.A[33]*params.A[33]+params.A[34]*params.A[34]+params.A[35]*params.A[35]+params.A[36]*params.A[36]+params.A[37]*params.A[37]+params.A[38]*params.A[38]+params.A[39]*params.A[39];
  work.quad_384992550912[0] = params.b[0]*params.b[0]+params.b[1]*params.b[1]+params.b[2]*params.b[2]+params.b[3]*params.b[3]+params.b[4]*params.b[4]+params.b[5]*params.b[5]+params.b[6]*params.b[6]+params.b[7]*params.b[7]+params.b[8]*params.b[8]+params.b[9]*params.b[9];
}
