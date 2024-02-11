/*
 * uavNavigationRL2D_2022a_icra.cpp
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "uavNavigationRL2D_2022a_icra".
 *
 * Model version              : 6.1
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C++ source code generated on : Tue Feb 28 16:30:11 2023
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "uavNavigationRL2D_2022a_icra.h"
#include "rtwtypes.h"
#include "uavNavigationRL2D_2022a_icra_private.h"
#include <cmath>
#include <emmintrin.h>
#include <cstring>
#include "rt_defines.h"

extern "C"
{

#include "rt_nonfinite.h"

}

/* Exported data definition */

/* Data with Exported storage */
real_T K_pitch{ 146.5804226 };         /* Referenced by: '<S5>/Gain10' */

real_T K_roll{ 150.1197777 };          /* Referenced by: '<S5>/Gain8' */

real_T K_z{ 0.839625266 };             /* Referenced by:
                                        * '<S5>/Constant'
                                        * '<S5>/Gain7'
                                        */

real_T T_aero{ 0.563387457 };          /* Referenced by:
                                        * '<S5>/Constant'
                                        * '<S5>/Gain7'
                                        */

real_T T_motor{ 0.049237556 };         /* Referenced by:
                                        * '<S3>/Gain'
                                        * '<S3>/Gain1'
                                        * '<S3>/Gain2'
                                        * '<S3>/Gain3'
                                        */

real_T kd_pitch{ 0.06518965284 };      /* Referenced by: '<S12>/Constant2' */

real_T kd_roll{ 0.06365268464 };       /* Referenced by: '<S13>/Constant4' */

real_T kd_x{ 2.605 };                  /* Referenced by: '<S15>/Constant5' */

real_T kd_y{ 2.725 };                  /* Referenced by: '<S16>/Constant4' */

real_T kd_yaw{ 0.12 };                 /* Referenced by: '<S14>/Gain2' */

real_T kd_z{ 8.14 };                   /* Referenced by: '<S17>/Constant4' */

real_T kp_pitch{ 0.4340612278 };       /* Referenced by: '<S12>/Constant1' */

real_T kp_roll{ 0.4238274211 };        /* Referenced by: '<S13>/Constant3' */

real_T kp_x{ 14.33 };                  /* Referenced by: '<S15>/Constant2' */

real_T kp_y{ 15.0 };                   /* Referenced by: '<S16>/Constant3' */

real_T kp_yaw{ 0.6 };                  /* Referenced by: '<S14>/Gain3' */

real_T kp_z{ 27.63 };                  /* Referenced by: '<S17>/Constant3' */

real_T max_u{ 1.0 };                   /* Referenced by:
                                        * '<S3>/Integrator Limited'
                                        * '<S3>/Integrator Limited1'
                                        * '<S3>/Integrator Limited2'
                                        * '<S3>/Integrator Limited3'
                                        */

real_T min_u{ -1.0 };                  /* Referenced by:
                                        * '<S3>/Integrator Limited'
                                        * '<S3>/Integrator Limited1'
                                        * '<S3>/Integrator Limited2'
                                        * '<S3>/Integrator Limited3'
                                        */

real_T tau_act{ 0.0 };                 /* Referenced by:
                                        * '<S3>/Transport Delay'
                                        * '<S3>/Transport Delay1'
                                        * '<S3>/Transport Delay2'
                                        * '<S3>/Transport Delay3'
                                        */

real_T tau_imu{ 0.05 };                /* Referenced by:
                                        * '<S5>/Transport Delay3'
                                        * '<S5>/Transport Delay4'
                                        * '<S5>/Transport Delay5'
                                        */

real_T tau_pitch{ 0.017482211 };       /* Referenced by:
                                        * '<S5>/Transport Delay11'
                                        * '<S5>/Transport Delay16'
                                        */

real_T tau_roll{ 0.017482211 };        /* Referenced by:
                                        * '<S5>/Transport Delay15'
                                        * '<S5>/Transport Delay6'
                                        */

real_T tau_x{ 0.1 };                 /* Referenced by: '<S5>/Transport Delay' */

real_T tau_y{ 0.1 };                /* Referenced by: '<S5>/Transport Delay1' */

real_T tau_z{ 0.05 };               /* Referenced by: '<S5>/Transport Delay2' */

real_T uavNavigati_velocity_references[3];/* '<Root>/velocity_references' */
real_T uavNavigationRL2D_2022a_u_cmd_z;/* '<Root>/u_cmd_z' */
real_T uavNavigationRL2D_states_output[12];/* '<Root>/states_output' */
real_T uavNavigationRL2_velocity_cmd_x;/* '<Root>/velocity_cmd_x' */
real_T uavNavigationRL2_velocity_cmd_y;/* '<Root>/velocity_cmd_y' */
real_T uavNavigationRL2_velocity_cmd_z;/* '<Root>/velocity_cmd_z' */
real_T uav_initial_x_position{ 0.0 };  /* Referenced by: '<S5>/Constant1' */

real_T uav_initial_y_position{ 0.0 };  /* Referenced by: '<S5>/Constant5' */

real_T uav_initial_z_position{ 0.0 };  /* Referenced by: '<S5>/Constant6' */

/* Forward declaration for local functions */
static void uavNavigationR_binary_expand_op(real_T in1[3], const int8_T
  in2_data[], const real_T in3_data[], const int32_T *in3_size, const real_T
  in4[4], const int32_T in5_size[2]);
static real_T uavNavigationRL2D_2022a_i_xnrm2(int32_T n, const real_T x[9],
  int32_T ix0);
static real_T uavNavigationRL2D_2022a_i_xdotc(int32_T n, const real_T x[9],
  int32_T ix0, const real_T y[9], int32_T iy0);
static void uavNavigationRL2D_2022a_i_xaxpy(int32_T n, real_T a, int32_T ix0,
  real_T y[9], int32_T iy0);
static real_T uavNavigationRL2D_2022a_xnrm2_p(const real_T x[3], int32_T ix0);
static void uavNavigationRL2D_2022a_i_xrotg(real_T *a, real_T *b, real_T *c,
  real_T *s);
static void uavNavigationRL2D_2022a_xaxpy_j(int32_T n, real_T a, const real_T x
  [9], int32_T ix0, real_T y[3], int32_T iy0);
static void uavNavigationRL2D_2022_xaxpy_jf(int32_T n, real_T a, const real_T x
  [3], int32_T ix0, real_T y[9], int32_T iy0);
static void uavNavigationRL2D_2022a_ic_xrot(real_T x[9], int32_T ix0, int32_T
  iy0, real_T c, real_T s);
static void uavNavigationRL2D_2022a_i_xswap(real_T x[9], int32_T ix0, int32_T
  iy0);
static void uavNavigationRL2D_2022a_icr_svd(const real_T A[9], real_T U[9],
  real_T s[3], real_T V[9]);
static real_T uavNavigationRL2D_2022a_ic_norm(const real_T x[2]);
static real_T uavNavigationRL2D_2022a__norm_h(const real_T x[3]);
static void rate_scheduler(RT_MODEL_uavNavigationRL2D_20_T *const
  uavNavigationRL2D_2022a_icra_M);

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 * (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 * (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                 /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;
  real_T* tBuf = uBuf + bufSz;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * std::abs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(RT_MODEL_uavNavigationRL2D_20_T *const
  uavNavigationRL2D_2022a_icra_M)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2])++;
  if ((uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2]) > 9) {/* Sample time: [0.005s, 0.0s] */
    uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE8 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si ,
  RT_MODEL_uavNavigationRL2D_20_T *const uavNavigationRL2D_2022a_icra_M)
{
  /* Solver Matrices */
#define uavNavigationRL2D_2022a_NSTAGES 13

  static real_T rt_ODE8_B[13]{
    4.174749114153025E-2, 0.0, 0.0, 0.0,
    0.0, -5.54523286112393E-2, 2.393128072011801E-1, 7.03510669403443E-1,
    -7.597596138144609E-1, 6.605630309222863E-1, 1.581874825101233E-1,
    -2.381095387528628E-1, 2.5E-1
  };

  static real_T rt_ODE8_C[13]{
    0.0, 5.555555555555556E-2, 8.333333333333333E-2, 1.25E-1,
    3.125E-1, 3.75E-1, 1.475E-1, 4.65E-1,
    5.648654513822596E-1, 6.5E-1, 9.246562776405044E-1, 1.0, 1.0
  };

  static real_T rt_ODE8_A[13][13]{
    /* rt_ODE8_A[0][] */
    { 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[1][] */
    { 5.555555555555556E-2, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[2][] */
    { 2.083333333333333E-2, 6.25E-2, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[3][] */
    { 3.125E-2, 0.0, 9.375E-2, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[4][] */
    { 3.125E-1, 0.0, -1.171875, 1.171875,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[5][] */
    { 3.75E-2, 0.0, 0.0, 1.875E-1,
      1.5E-1, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[6][] */
    { 4.791013711111111E-2, 0.0, 0.0, 1.122487127777778E-1,
      -2.550567377777778E-2, 1.284682388888889E-2, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[7][] */
    { 1.691798978729228E-2, 0.0, 0.0, 3.878482784860432E-1,
      3.597736985150033E-2, 1.969702142156661E-1, -1.727138523405018E-1, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[8][] */
    { 6.90957533591923E-2, 0.0, 0.0, -6.342479767288542E-1,
      -1.611975752246041E-1, 1.386503094588253E-1, 9.409286140357563E-1,
      2.11636326481944E-1,
      0.0, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[9][] */
    { 1.835569968390454E-1, 0.0, 0.0, -2.468768084315592,
      -2.912868878163005E-1, -2.647302023311738E-2, 2.8478387641928,
      2.813873314698498E-1,
      1.237448998633147E-1, 0.0, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[10][] */
    { -1.215424817395888, 0.0, 0.0, 1.667260866594577E1,
      9.157418284168179E-1, -6.056605804357471, -1.600357359415618E1,
      1.484930308629766E1,
      -1.337157573528985E1, 5.134182648179638, 0.0, 0.0, 0.0 },

    /* rt_ODE8_A[11][] */
    { 2.588609164382643E-1, 0.0, 0.0, -4.774485785489205,
      -4.350930137770325E-1, -3.049483332072241, 5.577920039936099,
      6.155831589861039,
      -5.062104586736938, 2.193926173180679, 1.346279986593349E-1, 0.0, 0.0 },

    /* rt_ODE8_A[12][] */
    { 8.224275996265075E-1, 0.0, 0.0, -1.165867325727766E1,
      -7.576221166909362E-1, 7.139735881595818E-1, 1.207577498689006E1,
      -2.127659113920403,
      1.990166207048956, -2.342864715440405E-1, 1.758985777079423E-1, 0.0, 0.0 },
  };

  time_T t { rtsiGetT(si) };

  time_T tnew { rtsiGetSolverStopTime(si) };

  time_T h { rtsiGetStepSize(si) };

  real_T *x { rtsiGetContStates(si) };

  ODE8_IntgData *intgData { static_cast<ODE8_IntgData *>(rtsiGetSolverData(si))
  };

  real_T *deltaY { intgData->deltaY };

  real_T *x0 { intgData->x0 };

  real_T* f[uavNavigationRL2D_2022a_NSTAGES];
  int idx,stagesIdx,statesIdx;
  double deltaX;
  int_T nXc { 14 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  f[0] = intgData->f[0];
  f[1] = intgData->f[1];
  f[2] = intgData->f[2];
  f[3] = intgData->f[3];
  f[4] = intgData->f[4];
  f[5] = intgData->f[5];
  f[6] = intgData->f[6];
  f[7] = intgData->f[7];
  f[8] = intgData->f[8];
  f[9] = intgData->f[9];
  f[10] = intgData->f[10];
  f[11] = intgData->f[11];
  f[12] = intgData->f[12];

  /* Save the state values at time t in y and x0*/
  (void) std::memset(deltaY, 0,
                     static_cast<uint_T>(nXc)*sizeof(real_T));
  (void) std::memcpy(x0, x,
                     nXc*sizeof(real_T));
  for (stagesIdx=0;stagesIdx<uavNavigationRL2D_2022a_NSTAGES;stagesIdx++) {
    (void) std::memcpy(x, x0,
                       static_cast<uint_T>(nXc)*sizeof(real_T));
    for (statesIdx=0;statesIdx<nXc;statesIdx++) {
      deltaX = 0;
      for (idx=0;idx<stagesIdx;idx++) {
        deltaX = deltaX + h*rt_ODE8_A[stagesIdx][idx]*f[idx][statesIdx];
      }

      x[statesIdx] = x0[statesIdx] + deltaX;
    }

    if (stagesIdx==0) {
      rtsiSetdX(si, f[stagesIdx]);
      uavNavigationRL2D_2022a_icra_derivatives(uavNavigationRL2D_2022a_icra_M);
    } else {
      (stagesIdx==uavNavigationRL2D_2022a_NSTAGES-1)? rtsiSetT(si, tnew) :
        rtsiSetT(si, t + h*rt_ODE8_C[stagesIdx]);
      rtsiSetdX(si, f[stagesIdx]);
      uavNavigationRL2D_2022a_icra_step(uavNavigationRL2D_2022a_icra_M);
      uavNavigationRL2D_2022a_icra_derivatives(uavNavigationRL2D_2022a_icra_M);
    }

    for (statesIdx=0;statesIdx<nXc;statesIdx++) {
      deltaY[statesIdx] = deltaY[statesIdx] + h*rt_ODE8_B[stagesIdx]*f[stagesIdx]
        [statesIdx];
    }
  }

  for (statesIdx=0;statesIdx<nXc;statesIdx++) {
    x[statesIdx] = x0[statesIdx] + deltaY[statesIdx];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S19>/reference correction'
 *    '<S20>/reference correction'
 *    '<S21>/reference correction'
 *    '<S34>/reference correction'
 */
void uavNavigati_referencecorrection(real_T rtu_yaw, real_T rtu_x, real_T rtu_y,
  B_referencecorrection_uavNavi_T *localB)
{
  real_T out_idx_0;
  real_T out_idx_1;
  real_T tmp;
  real_T tmp_0;
  real_T tmp_1;
  out_idx_0 = std::cos(rtu_yaw);
  tmp = std::sin(rtu_yaw);
  tmp_0 = std::sin(rtu_yaw);
  tmp_1 = std::cos(rtu_yaw);
  tmp_0 = -tmp_0;
  out_idx_1 = out_idx_0 * rtu_x;
  out_idx_1 += tmp * rtu_y;
  out_idx_0 = out_idx_1;
  out_idx_1 = tmp_0 * rtu_x;
  out_idx_1 += tmp_1 * rtu_y;
  localB->x_out = out_idx_0;
  localB->y_out = out_idx_1;
}

/*
 * Termination for atomic system:
 *    '<S19>/reference correction'
 *    '<S20>/reference correction'
 *    '<S21>/reference correction'
 *    '<S34>/reference correction'
 */
void uavNav_referencecorrection_Term(void)
{
}

void rt_mrdivided4x4_snf(const real_T u0[16], const real_T u1[16], real_T y[16])
{
  real_T A[16];
  real_T s;
  real_T smax;
  int32_T c;
  int32_T ia0;
  int32_T jBcol;
  int32_T jj;
  int32_T jm1;
  int32_T jpiv;
  int32_T mmj;
  int8_T ipiv[4];
  std::memcpy(&y[0], &u0[0], sizeof(real_T) << 4U);
  std::memcpy(&A[0], &u1[0], sizeof(real_T) << 4U);
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  for (int32_T ib0{0}; ib0 < 3; ib0++) {
    real_T x;
    int32_T ONE;
    int32_T ix;
    int32_T jp1j;
    ia0 = ib0 + 1;
    jm1 = ia0 - 1;
    mmj = 4 - ia0;
    c = jm1 * 5;
    jj = c + 1;
    jp1j = jj + 1;
    c = mmj + 1;
    ONE = 1;
    ix = jj - 1;
    x = A[jj - 1];
    s = std::abs(x);
    smax = s;
    for (jpiv = 2; jpiv <= c; jpiv++) {
      ix++;
      x = A[ix];
      s = std::abs(x);
      if (s > smax) {
        ONE = jpiv;
        smax = s;
      }
    }

    ONE--;
    jpiv = (jj + ONE) - 1;
    if (A[jpiv] != 0.0) {
      if (ONE != 0) {
        c = ia0 + ONE;
        ipiv[ia0 - 1] = static_cast<int8_T>(c);
        ONE += jm1;
        smax = A[jm1];
        A[jm1] = A[ONE];
        A[ONE] = smax;
        jm1 += 4;
        ONE += 4;
        smax = A[jm1];
        A[jm1] = A[ONE];
        A[ONE] = smax;
        jm1 += 4;
        ONE += 4;
        smax = A[jm1];
        A[jm1] = A[ONE];
        A[ONE] = smax;
        jm1 += 4;
        ONE += 4;
        smax = A[jm1];
        A[jm1] = A[ONE];
        A[ONE] = smax;
      }

      c = mmj - 1;
      jm1 = jp1j + c;
      for (c = jp1j; c <= jm1; c++) {
        x = A[c - 1];
        s = A[jj - 1];
        s = x / s;
        A[c - 1] = s;
      }
    }

    c = 3 - ia0;
    ONE = jj + 3;
    ia0 = jj + 5;
    jpiv = ia0 - 1;
    for (ia0 = 0; ia0 <= c; ia0++) {
      smax = A[ONE];
      if (smax != 0.0) {
        smax = -smax;
        ix = jp1j - 1;
        jm1 = jpiv + 1;
        jj = mmj + jpiv;
        for (jBcol = jm1; jBcol <= jj; jBcol++) {
          A[jBcol - 1] += A[ix] * smax;
          ix++;
        }
      }

      ONE += 4;
      jpiv += 4;
    }
  }

  for (int32_T ib0{0}; ib0 < 4; ib0++) {
    ia0 = ib0 + 1;
    jBcol = ((ia0 - 1) << 2) - 1;
    mmj = ((ia0 - 1) << 2) - 1;
    jm1 = ia0 - 2;
    for (jpiv = 0; jpiv <= jm1; jpiv++) {
      c = jpiv + 1;
      jj = ((c - 1) << 2) - 1;
      if (A[c + mmj] != 0.0) {
        y[jBcol + 1] -= A[c + mmj] * y[jj + 1];
        y[jBcol + 2] -= A[c + mmj] * y[jj + 2];
        y[jBcol + 3] -= A[c + mmj] * y[jj + 3];
        y[jBcol + 4] -= A[c + mmj] * y[jj + 4];
      }
    }

    s = A[ia0 + mmj];
    smax = 1.0 / s;
    y[jBcol + 1] *= smax;
    y[jBcol + 2] *= smax;
    y[jBcol + 3] *= smax;
    y[jBcol + 4] *= smax;
  }

  for (int32_T ib0{3}; ib0 >= 0; ib0--) {
    jBcol = (ib0 << 2) - 1;
    mmj = (ib0 << 2) - 1;
    jm1 = ib0 + 2;
    for (jpiv = jm1; jpiv < 5; jpiv++) {
      jj = ((jpiv - 1) << 2) - 1;
      if (A[jpiv + mmj] != 0.0) {
        y[jBcol + 1] -= A[jpiv + mmj] * y[jj + 1];
        y[jBcol + 2] -= A[jpiv + mmj] * y[jj + 2];
        y[jBcol + 3] -= A[jpiv + mmj] * y[jj + 3];
        y[jBcol + 4] -= A[jpiv + mmj] * y[jj + 4];
      }
    }
  }

  for (int32_T ib0{2}; ib0 >= 0; ib0--) {
    int8_T ipiv_0;
    ipiv_0 = ipiv[ib0];
    if (ib0 + 1 != ipiv_0) {
      ia0 = ipiv_0 - 1;
      smax = y[ib0 << 2];
      y[ib0 << 2] = y[ia0 << 2];
      y[ia0 << 2] = smax;
      smax = y[(ib0 << 2) + 1];
      y[(ib0 << 2) + 1] = y[(ia0 << 2) + 1];
      y[(ia0 << 2) + 1] = smax;
      smax = y[(ib0 << 2) + 2];
      y[(ib0 << 2) + 2] = y[(ia0 << 2) + 2];
      y[(ia0 << 2) + 2] = smax;
      smax = y[(ib0 << 2) + 3];
      y[(ib0 << 2) + 3] = y[(ia0 << 2) + 3];
      y[(ia0 << 2) + 3] = smax;
    }
  }
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int32_T tmp;
    int32_T tmp_0;
    if (u1 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u0 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = std::atan2(static_cast<real_T>(tmp_0), static_cast<real_T>(tmp));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }

  return y;
}

static void uavNavigationR_binary_expand_op(real_T in1[3], const int8_T
  in2_data[], const real_T in3_data[], const int32_T *in3_size, const real_T
  in4[4], const int32_T in5_size[2])
{
  int32_T in5_idx_0;
  int32_T loop_ub;

  /* MATLAB Function: '<S41>/kinematic differential equation' */
  in5_idx_0 = in5_size[1];
  loop_ub = in5_idx_0 == 1 ? *in3_size : in5_idx_0;
  for (in5_idx_0 = 0; in5_idx_0 < loop_ub; in5_idx_0++) {
    real_T varargin_1;
    real_T varargin_2;
    varargin_1 = in4[1];
    varargin_2 = in4[0];
    varargin_1 = rt_atan2d_snf(varargin_1, varargin_2);
    in1[in2_data[0] - 1] = -in3_data[0] * 2.0 * varargin_1;
  }

  /* End of MATLAB Function: '<S41>/kinematic differential equation' */
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static real_T uavNavigationRL2D_2022a_i_xnrm2(int32_T n, const real_T x[9],
  int32_T ix0)
{
  real_T scale;
  real_T y;
  int32_T kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = (ix0 + n) - 1;
  for (int32_T k{ix0}; k <= kend; k++) {
    real_T absxk;
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      real_T t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      real_T t;
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static real_T uavNavigationRL2D_2022a_i_xdotc(int32_T n, const real_T x[9],
  int32_T ix0, const real_T y[9], int32_T iy0)
{
  real_T d;
  int32_T b;
  int32_T ix;
  int32_T iy;
  ix = ix0 - 1;
  iy = iy0 - 1;
  d = 0.0;
  b = static_cast<uint8_T>(n);
  for (int32_T k{0}; k < b; k++) {
    d += x[ix + k] * y[iy + k];
  }

  return d;
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_i_xaxpy(int32_T n, real_T a, int32_T ix0,
  real_T y[9], int32_T iy0)
{
  if (!(a == 0.0)) {
    int32_T iy;
    iy = iy0 - 1;
    for (int32_T k{0}; k < n; k++) {
      y[iy + k] += y[(ix0 + k) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static real_T uavNavigationRL2D_2022a_xnrm2_p(const real_T x[3], int32_T ix0)
{
  real_T scale;
  real_T y;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  for (int32_T k{ix0}; k <= ix0 + 1; k++) {
    real_T absxk;
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      real_T t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      real_T t;
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_i_xrotg(real_T *a, real_T *b, real_T *c,
  real_T *s)
{
  real_T absa;
  real_T absb;
  real_T roe;
  real_T scale;
  roe = *b;
  absa = std::abs(*a);
  absb = std::abs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    real_T ads;
    real_T bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = scale;
  }
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_xaxpy_j(int32_T n, real_T a, const real_T x
  [9], int32_T ix0, real_T y[3], int32_T iy0)
{
  if (!(a == 0.0)) {
    int32_T iy;
    int32_T scalarLB;
    int32_T vectorUB;
    iy = iy0 - 1;
    scalarLB = (n / 2) << 1;
    vectorUB = scalarLB - 2;
    for (int32_T k{0}; k <= vectorUB; k += 2) {
      __m128d tmp;
      __m128d tmp_0;
      tmp = _mm_loadu_pd(&x[(ix0 + k) - 1]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(a));
      tmp_0 = _mm_loadu_pd(&y[iy + k]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[iy + k], tmp);
    }

    for (int32_T k{scalarLB}; k < n; k++) {
      y[iy + k] += x[(ix0 + k) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022_xaxpy_jf(int32_T n, real_T a, const real_T x
  [3], int32_T ix0, real_T y[9], int32_T iy0)
{
  if (!(a == 0.0)) {
    int32_T iy;
    int32_T scalarLB;
    int32_T vectorUB;
    iy = iy0 - 1;
    scalarLB = (n / 2) << 1;
    vectorUB = scalarLB - 2;
    for (int32_T k{0}; k <= vectorUB; k += 2) {
      __m128d tmp;
      __m128d tmp_0;
      tmp = _mm_loadu_pd(&x[(ix0 + k) - 1]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(a));
      tmp_0 = _mm_loadu_pd(&y[iy + k]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&y[iy + k], tmp);
    }

    for (int32_T k{scalarLB}; k < n; k++) {
      y[iy + k] += x[(ix0 + k) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_ic_xrot(real_T x[9], int32_T ix0, int32_T
  iy0, real_T c, real_T s)
{
  real_T temp;
  int32_T ix;
  int32_T iy;
  ix = ix0 - 1;
  iy = iy0 - 1;
  temp = c * x[ix] + s * x[iy];
  x[iy] = c * x[iy] - s * x[ix];
  x[ix] = temp;
  temp = x[ix + 1] * c + x[iy + 1] * s;
  x[iy + 1] = x[iy + 1] * c - x[ix + 1] * s;
  x[ix + 1] = temp;
  temp = x[ix + 2] * c + x[iy + 2] * s;
  x[iy + 2] = x[iy + 2] * c - x[ix + 2] * s;
  x[ix + 2] = temp;
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_i_xswap(real_T x[9], int32_T ix0, int32_T
  iy0)
{
  real_T temp;
  int32_T ix;
  int32_T iy;
  ix = ix0 - 1;
  iy = iy0 - 1;
  temp = x[ix];
  x[ix] = x[iy];
  x[iy] = temp;
  temp = x[ix + 1];
  x[ix + 1] = x[iy + 1];
  x[iy + 1] = temp;
  temp = x[ix + 2];
  x[ix + 2] = x[iy + 2];
  x[iy + 2] = temp;
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static void uavNavigationRL2D_2022a_icr_svd(const real_T A[9], real_T U[9],
  real_T s[3], real_T V[9])
{
  __m128d tmp;
  real_T b_A[9];
  real_T b_s[3];
  real_T e[3];
  real_T work[3];
  real_T nrm;
  real_T rt;
  real_T smm1;
  real_T sqds;
  real_T ztest;
  int32_T kase;
  int32_T m;
  int32_T qjj;
  int32_T qp1;
  int32_T qq;
  int32_T scalarLB;
  int32_T vectorUB;
  b_s[0] = 0.0;
  e[0] = 0.0;
  work[0] = 0.0;
  b_s[1] = 0.0;
  e[1] = 0.0;
  work[1] = 0.0;
  b_s[2] = 0.0;
  e[2] = 0.0;
  work[2] = 0.0;
  for (m = 0; m < 9; m++) {
    b_A[m] = A[m];
    U[m] = 0.0;
    V[m] = 0.0;
  }

  for (m = 0; m < 2; m++) {
    boolean_T apply_transform;
    qp1 = m + 2;
    qq = (3 * m + m) + 1;
    apply_transform = false;
    nrm = uavNavigationRL2D_2022a_i_xnrm2(3 - m, b_A, qq);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qq - 1] < 0.0) {
        b_s[m] = -nrm;
      } else {
        b_s[m] = nrm;
      }

      if (std::abs(b_s[m]) >= 1.0020841800044864E-292) {
        nrm = 1.0 / b_s[m];
        qjj = (qq - m) + 2;
        scalarLB = ((((qjj - qq) + 1) / 2) << 1) + qq;
        vectorUB = scalarLB - 2;
        for (kase = qq; kase <= vectorUB; kase += 2) {
          tmp = _mm_loadu_pd(&b_A[kase - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(nrm));
          _mm_storeu_pd(&b_A[kase - 1], tmp);
        }

        for (kase = scalarLB; kase <= qjj; kase++) {
          b_A[kase - 1] *= nrm;
        }
      } else {
        qjj = (qq - m) + 2;
        scalarLB = ((((qjj - qq) + 1) / 2) << 1) + qq;
        vectorUB = scalarLB - 2;
        for (kase = qq; kase <= vectorUB; kase += 2) {
          tmp = _mm_loadu_pd(&b_A[kase - 1]);
          tmp = _mm_div_pd(tmp, _mm_set1_pd(b_s[m]));
          _mm_storeu_pd(&b_A[kase - 1], tmp);
        }

        for (kase = scalarLB; kase <= qjj; kase++) {
          b_A[kase - 1] /= b_s[m];
        }
      }

      b_A[qq - 1]++;
      b_s[m] = -b_s[m];
    } else {
      b_s[m] = 0.0;
    }

    for (kase = qp1; kase < 4; kase++) {
      qjj = (kase - 1) * 3 + m;
      if (apply_transform) {
        uavNavigationRL2D_2022a_i_xaxpy(3 - m, -(uavNavigationRL2D_2022a_i_xdotc
          (3 - m, b_A, qq, b_A, qjj + 1) / b_A[m + 3 * m]), qq, b_A, qjj + 1);
      }

      e[kase - 1] = b_A[qjj];
    }

    for (qq = m + 1; qq < 4; qq++) {
      U[(qq + 3 * m) - 1] = b_A[(3 * m + qq) - 1];
    }

    if (m + 1 <= 1) {
      nrm = uavNavigationRL2D_2022a_xnrm2_p(e, 2);
      if (nrm == 0.0) {
        e[0] = 0.0;
      } else {
        if (e[1] < 0.0) {
          e[0] = -nrm;
        } else {
          e[0] = nrm;
        }

        nrm = e[0];
        if (std::abs(e[0]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[0];
          scalarLB = (((4 - qp1) / 2) << 1) + qp1;
          vectorUB = scalarLB - 2;
          for (qq = qp1; qq <= vectorUB; qq += 2) {
            tmp = _mm_loadu_pd(&e[qq - 1]);
            tmp = _mm_mul_pd(tmp, _mm_set1_pd(nrm));
            _mm_storeu_pd(&e[qq - 1], tmp);
          }

          for (qq = scalarLB; qq < 4; qq++) {
            e[qq - 1] *= nrm;
          }
        } else {
          scalarLB = (((4 - qp1) / 2) << 1) + qp1;
          vectorUB = scalarLB - 2;
          for (qq = qp1; qq <= vectorUB; qq += 2) {
            tmp = _mm_loadu_pd(&e[qq - 1]);
            tmp = _mm_div_pd(tmp, _mm_set1_pd(nrm));
            _mm_storeu_pd(&e[qq - 1], tmp);
          }

          for (qq = scalarLB; qq < 4; qq++) {
            e[qq - 1] /= nrm;
          }
        }

        e[1]++;
        e[0] = -e[0];
        for (qq = qp1; qq < 4; qq++) {
          work[qq - 1] = 0.0;
        }

        for (qq = qp1; qq < 4; qq++) {
          uavNavigationRL2D_2022a_xaxpy_j(2, e[qq - 1], b_A, 3 * (qq - 1) + 2,
            work, 2);
        }

        for (qq = qp1; qq < 4; qq++) {
          uavNavigationRL2D_2022_xaxpy_jf(2, -e[qq - 1] / e[1], work, 2, b_A, 3 *
            (qq - 1) + 2);
        }
      }

      for (qq = qp1; qq < 4; qq++) {
        V[qq - 1] = e[qq - 1];
      }
    }
  }

  m = 1;
  b_s[2] = b_A[8];
  e[1] = b_A[7];
  e[2] = 0.0;
  U[6] = 0.0;
  U[7] = 0.0;
  U[8] = 1.0;
  for (qp1 = 1; qp1 >= 0; qp1--) {
    qq = 3 * qp1 + qp1;
    if (b_s[qp1] != 0.0) {
      for (kase = qp1 + 2; kase < 4; kase++) {
        qjj = ((kase - 1) * 3 + qp1) + 1;
        uavNavigationRL2D_2022a_i_xaxpy(3 - qp1,
          -(uavNavigationRL2D_2022a_i_xdotc(3 - qp1, U, qq + 1, U, qjj) / U[qq]),
          qq + 1, U, qjj);
      }

      for (qjj = qp1 + 1; qjj < 4; qjj++) {
        U[(qjj + 3 * qp1) - 1] = -U[(3 * qp1 + qjj) - 1];
      }

      U[qq]++;
      qq = qp1;
      if (qq - 1 >= 0) {
        U[3 * qp1] = 0.0;
      }
    } else {
      U[3 * qp1] = 0.0;
      U[3 * qp1 + 1] = 0.0;
      U[3 * qp1 + 2] = 0.0;
      U[qq] = 1.0;
    }
  }

  for (qp1 = 2; qp1 >= 0; qp1--) {
    if ((qp1 + 1 <= 1) && (e[0] != 0.0)) {
      uavNavigationRL2D_2022a_i_xaxpy(2, -(uavNavigationRL2D_2022a_i_xdotc(2, V,
        2, V, 5) / V[1]), 2, V, 5);
      uavNavigationRL2D_2022a_i_xaxpy(2, -(uavNavigationRL2D_2022a_i_xdotc(2, V,
        2, V, 8) / V[1]), 2, V, 8);
    }

    V[3 * qp1] = 0.0;
    V[3 * qp1 + 1] = 0.0;
    V[3 * qp1 + 2] = 0.0;
    V[qp1 + 3 * qp1] = 1.0;
  }

  for (qp1 = 0; qp1 < 3; qp1++) {
    smm1 = e[qp1];
    if (b_s[qp1] != 0.0) {
      rt = std::abs(b_s[qp1]);
      nrm = b_s[qp1] / rt;
      b_s[qp1] = rt;
      if (qp1 + 1 < 3) {
        smm1 /= nrm;
      }

      qq = 3 * qp1 + 1;
      scalarLB = 2 + qq;
      vectorUB = scalarLB - 2;
      for (qjj = qq; qjj <= vectorUB; qjj += 2) {
        tmp = _mm_loadu_pd(&U[qjj - 1]);
        tmp = _mm_mul_pd(tmp, _mm_set1_pd(nrm));
        _mm_storeu_pd(&U[qjj - 1], tmp);
      }

      for (qjj = scalarLB; qjj <= qq + 2; qjj++) {
        U[qjj - 1] *= nrm;
      }
    }

    if ((qp1 + 1 < 3) && (smm1 != 0.0)) {
      rt = std::abs(smm1);
      nrm = rt / smm1;
      smm1 = rt;
      b_s[qp1 + 1] *= nrm;
      qq = (qp1 + 1) * 3 + 1;
      scalarLB = 2 + qq;
      vectorUB = scalarLB - 2;
      for (qjj = qq; qjj <= vectorUB; qjj += 2) {
        tmp = _mm_loadu_pd(&V[qjj - 1]);
        tmp = _mm_mul_pd(tmp, _mm_set1_pd(nrm));
        _mm_storeu_pd(&V[qjj - 1], tmp);
      }

      for (qjj = scalarLB; qjj <= qq + 2; qjj++) {
        V[qjj - 1] *= nrm;
      }
    }

    e[qp1] = smm1;
  }

  qp1 = 0;
  nrm = std::fmax(0.0, std::fmax(std::abs(b_s[0]), std::abs(e[0])));
  nrm = std::fmax(nrm, std::fmax(std::abs(b_s[1]), std::abs(e[1])));
  nrm = std::fmax(nrm, std::fmax(std::abs(b_s[2]), std::abs(e[2])));
  while ((m + 2 > 0) && (qp1 < 75)) {
    kase = m + 1;
    int32_T exitg1;
    do {
      exitg1 = 0;
      qq = kase;
      if (kase == 0) {
        exitg1 = 1;
      } else {
        rt = std::abs(e[kase - 1]);
        if ((rt <= (std::abs(b_s[kase - 1]) + std::abs(b_s[kase])) *
             2.2204460492503131E-16) || ((rt <= 1.0020841800044864E-292) ||
             ((qp1 > 20) && (rt <= 2.2204460492503131E-16 * nrm)))) {
          e[kase - 1] = 0.0;
          exitg1 = 1;
        } else {
          kase--;
        }
      }
    } while (exitg1 == 0);

    if (m + 1 == kase) {
      kase = 4;
    } else {
      boolean_T exitg2;
      qjj = m + 2;
      scalarLB = m + 2;
      exitg2 = false;
      while ((!exitg2) && (scalarLB >= kase)) {
        qjj = scalarLB;
        if (scalarLB == kase) {
          exitg2 = true;
        } else {
          rt = 0.0;
          if (scalarLB < m + 2) {
            rt = std::abs(e[scalarLB - 1]);
          }

          if (scalarLB > kase + 1) {
            rt += std::abs(e[scalarLB - 2]);
          }

          ztest = std::abs(b_s[scalarLB - 1]);
          if ((ztest <= 2.2204460492503131E-16 * rt) || (ztest <=
               1.0020841800044864E-292)) {
            b_s[scalarLB - 1] = 0.0;
            exitg2 = true;
          } else {
            scalarLB--;
          }
        }
      }

      if (qjj == kase) {
        kase = 3;
      } else if (m + 2 == qjj) {
        kase = 1;
      } else {
        kase = 2;
        qq = qjj;
      }
    }

    switch (kase) {
     case 1:
      rt = e[m];
      e[m] = 0.0;
      for (qjj = m + 1; qjj >= qq + 1; qjj--) {
        smm1 = e[0];
        uavNavigationRL2D_2022a_i_xrotg(&b_s[qjj - 1], &rt, &ztest, &sqds);
        if (qjj > qq + 1) {
          rt = -sqds * smm1;
          smm1 *= ztest;
        }

        uavNavigationRL2D_2022a_ic_xrot(V, 3 * (qjj - 1) + 1, 3 * (m + 1) + 1,
          ztest, sqds);
        e[0] = smm1;
      }
      break;

     case 2:
      rt = e[qq - 1];
      e[qq - 1] = 0.0;
      for (qjj = qq + 1; qjj <= m + 2; qjj++) {
        uavNavigationRL2D_2022a_i_xrotg(&b_s[qjj - 1], &rt, &ztest, &sqds);
        rt = e[qjj - 1] * -sqds;
        e[qjj - 1] *= ztest;
        uavNavigationRL2D_2022a_ic_xrot(U, 3 * (qjj - 1) + 1, 3 * (qq - 1) + 1,
          ztest, sqds);
      }
      break;

     case 3:
      {
        real_T emm1;
        real_T shift;
        ztest = std::fmax(std::fmax(std::fmax(std::fmax(std::abs(b_s[m + 1]),
          std::abs(b_s[m])), std::abs(e[m])), std::abs(b_s[qq])), std::abs(e[qq]));
        rt = b_s[m + 1] / ztest;
        smm1 = b_s[m] / ztest;
        emm1 = e[m] / ztest;
        sqds = b_s[qq] / ztest;
        smm1 = ((smm1 + rt) * (smm1 - rt) + emm1 * emm1) / 2.0;
        emm1 *= rt;
        emm1 *= emm1;
        if ((smm1 != 0.0) || (emm1 != 0.0)) {
          shift = std::sqrt(smm1 * smm1 + emm1);
          if (smm1 < 0.0) {
            shift = -shift;
          }

          shift = emm1 / (smm1 + shift);
        } else {
          shift = 0.0;
        }

        rt = (sqds + rt) * (sqds - rt) + shift;
        ztest = e[qq] / ztest * sqds;
        for (qjj = qq + 1; qjj <= m + 1; qjj++) {
          uavNavigationRL2D_2022a_i_xrotg(&rt, &ztest, &sqds, &smm1);
          if (qjj > qq + 1) {
            e[0] = rt;
          }

          rt = b_s[qjj - 1] * sqds + e[qjj - 1] * smm1;
          e[qjj - 1] = e[qjj - 1] * sqds - b_s[qjj - 1] * smm1;
          ztest = smm1 * b_s[qjj];
          b_s[qjj] *= sqds;
          uavNavigationRL2D_2022a_ic_xrot(V, 3 * (qjj - 1) + 1, 3 * qjj + 1,
            sqds, smm1);
          b_s[qjj - 1] = rt;
          uavNavigationRL2D_2022a_i_xrotg(&b_s[qjj - 1], &ztest, &sqds, &smm1);
          rt = e[qjj - 1] * sqds + smm1 * b_s[qjj];
          b_s[qjj] = e[qjj - 1] * -smm1 + sqds * b_s[qjj];
          ztest = smm1 * e[qjj];
          e[qjj] *= sqds;
          uavNavigationRL2D_2022a_ic_xrot(U, 3 * (qjj - 1) + 1, 3 * qjj + 1,
            sqds, smm1);
        }

        e[m] = rt;
        qp1++;
      }
      break;

     default:
      if (b_s[qq] < 0.0) {
        b_s[qq] = -b_s[qq];
        qp1 = 3 * qq + 1;
        scalarLB = 2 + qp1;
        vectorUB = scalarLB - 2;
        for (qjj = qp1; qjj <= vectorUB; qjj += 2) {
          tmp = _mm_loadu_pd(&V[qjj - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(-1.0));
          _mm_storeu_pd(&V[qjj - 1], tmp);
        }

        for (qjj = scalarLB; qjj <= qp1 + 2; qjj++) {
          V[qjj - 1] = -V[qjj - 1];
        }
      }

      qp1 = qq + 1;
      while ((qq + 1 < 3) && (b_s[qq] < b_s[qp1])) {
        rt = b_s[qq];
        b_s[qq] = b_s[qp1];
        b_s[qp1] = rt;
        uavNavigationRL2D_2022a_i_xswap(V, 3 * qq + 1, 3 * (qq + 1) + 1);
        uavNavigationRL2D_2022a_i_xswap(U, 3 * qq + 1, 3 * (qq + 1) + 1);
        qq = qp1;
        qp1++;
      }

      qp1 = 0;
      m--;
      break;
    }
  }

  s[0] = b_s[0];
  s[1] = b_s[1];
  s[2] = b_s[2];
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static real_T uavNavigationRL2D_2022a_ic_norm(const real_T x[2])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S33>/attitude mapper ' */
static real_T uavNavigationRL2D_2022a__norm_h(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

/* Model step function */
void uavNavigationRL2D_2022a_icra_step(RT_MODEL_uavNavigationRL2D_20_T *const
  uavNavigationRL2D_2022a_icra_M)
{
  B_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_B{
    uavNavigationRL2D_2022a_icra_M->blockIO };

  DW_uavNavigationRL2D_2022a_ic_T *uavNavigationRL2D_2022a_icra_DW{
    uavNavigationRL2D_2022a_icra_M->dwork };

  X_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_X{
    uavNavigationRL2D_2022a_icra_M->contStates };

  if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M)) {
    /* set solver stop time */
    if (!(uavNavigationRL2D_2022a_icra_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(uavNavigationRL2D_2022a_icra_M->solverInfo,
                            ((uavNavigationRL2D_2022a_icra_M->Timing.clockTickH0
        + 1) * uavNavigationRL2D_2022a_icra_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(uavNavigationRL2D_2022a_icra_M->solverInfo,
                            ((uavNavigationRL2D_2022a_icra_M->Timing.clockTick0
        + 1) * uavNavigationRL2D_2022a_icra_M->Timing.stepSize0 +
        uavNavigationRL2D_2022a_icra_M->Timing.clockTickH0 *
        uavNavigationRL2D_2022a_icra_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(uavNavigationRL2D_2022a_icra_M)) {
    uavNavigationRL2D_2022a_icra_M->Timing.t[0] = rtsiGetT
      (uavNavigationRL2D_2022a_icra_M->solverInfo);
  }

  {
    static const int8_T c[3]{ 0, 0, 1 };

    __m128d tmp_0;
    __m128d tmp_1;
    real_T tmp[12];
    real_T a[9];
    real_T d_b[9];
    real_T d_y[9];
    real_T roll_pitch_DCM[9];
    real_T tempR[9];
    real_T b_y[4];
    real_T ct[3];
    real_T euler_angle[3];
    real_T out[2];
    const real_T *tmp_2;
    real_T absxk;
    real_T q_idx_1;
    real_T q_idx_2;
    real_T q_idx_3;
    real_T scale;
    real_T st_idx_1;
    real_T st_idx_2;
    real_T t;
    real_T x_data;
    real_T *lastU;
    int32_T f_size[2];
    int32_T i;
    int32_T x_size;
    int_T iy;
    int8_T g_data;
    boolean_T exitg1;
    boolean_T mask1;

    /* Integrator: '<S5>/Integrator7' */
    uavNavigationRL2D_2022a_icra_B->Integrator7 =
      uavNavigationRL2D_2022a_icra_X->Integrator7_CSTATE;

    /* Integrator: '<S5>/Integrator6' */
    uavNavigationRL2D_2022a_icra_B->Integrator6 =
      uavNavigationRL2D_2022a_icra_X->Integrator6_CSTATE;

    /* Integrator: '<S5>/Integrator5' */
    uavNavigationRL2D_2022a_icra_B->Integrator5 =
      uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE;

    /* Sum: '<S5>/Sum1' incorporates:
     *  Constant: '<S5>/Constant1'
     *  Constant: '<S5>/Constant5'
     *  Constant: '<S5>/Constant6'
     */
    uavNavigationRL2D_2022a_icra_B->Sum1[0] =
      uavNavigationRL2D_2022a_icra_B->Integrator7 + uav_initial_x_position;
    uavNavigationRL2D_2022a_icra_B->Sum1[1] =
      uavNavigationRL2D_2022a_icra_B->Integrator6 + uav_initial_y_position;
    uavNavigationRL2D_2022a_icra_B->Sum1[2] =
      uavNavigationRL2D_2022a_icra_B->Integrator5 + uav_initial_z_position;

    /* TransportDelay: '<S5>/Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_x;
      uavNavigationRL2D_2022a_icra_B->TransportDelay = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_y;
      uavNavigationRL2D_2022a_icra_B->TransportDelay1 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay1_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay2' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_z;
      uavNavigationRL2D_2022a_icra_B->TransportDelay2 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay2_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay3' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_imu;
      uavNavigationRL2D_2022a_icra_B->TransportDelay3 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay3_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay4' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay4_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_imu;
      uavNavigationRL2D_2022a_icra_B->TransportDelay4 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay4_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay5' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay5_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_imu;
      uavNavigationRL2D_2022a_icra_B->TransportDelay5 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay5_InitOutput,
        0,
        0);
    }

    /* TransportDelay: '<S5>/Transport Delay6' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay6_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_roll;
      uavNavigationRL2D_2022a_icra_B->TransportDelay6 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay6_InitOutput,
        1,
        0);
    }

    /* Gain: '<S6>/Gain' */
    uavNavigationRL2D_2022a_icra_B->Gain =
      uavNavigationRL2D_2022a_icra_P.Gain_Gain *
      uavNavigationRL2D_2022a_icra_B->TransportDelay6;

    /* Gain: '<S2>/Gain' */
    uavNavigationRL2D_2022a_icra_B->roll =
      uavNavigationRL2D_2022a_icra_P.Gain_Gain_h *
      uavNavigationRL2D_2022a_icra_B->Gain;

    /* TransportDelay: '<S5>/Transport Delay11' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay11_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_pitch;
      uavNavigationRL2D_2022a_icra_B->TransportDelay11 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay11_InitOutput,
        1,
        0);
    }

    /* Gain: '<S7>/Gain' */
    uavNavigationRL2D_2022a_icra_B->Gain_n =
      uavNavigationRL2D_2022a_icra_P.Gain_Gain_j *
      uavNavigationRL2D_2022a_icra_B->TransportDelay11;

    /* TransferFcn: '<S5>/Transfer Fcn7' */
    uavNavigationRL2D_2022a_icra_B->TransferFcn7 = 0.0;
    uavNavigationRL2D_2022a_icra_B->TransferFcn7 +=
      uavNavigationRL2D_2022a_icra_P.TransferFcn7_C *
      uavNavigationRL2D_2022a_icra_X->TransferFcn7_CSTATE;

    /* TransferFcn: '<S5>/Transfer Fcn8' */
    uavNavigationRL2D_2022a_icra_B->TransferFcn8 = 0.0;
    uavNavigationRL2D_2022a_icra_B->TransferFcn8 +=
      uavNavigationRL2D_2022a_icra_P.TransferFcn8_C *
      uavNavigationRL2D_2022a_icra_X->TransferFcn8_CSTATE;
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S5>/Matrix Divide' incorporates:
       *  Constant: '<S5>/Constant2'
       *  Constant: '<S5>/Constant4'
       */
      rt_mrdivided4x4_snf(uavNavigationRL2D_2022a_icra_P.Constant2_Value,
                          uavNavigationRL2D_2022a_icra_P.Constant4_Value,
                          uavNavigationRL2D_2022a_icra_B->MatrixDivide);
    }

    /* Integrator: '<S3>/Integrator Limited' */
    /* Limited  Integrator  */
    if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE >= max_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE = max_u;
    } else if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE <= min_u)
    {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE = min_u;
    }

    /* Integrator: '<S3>/Integrator Limited' */
    uavNavigationRL2D_2022a_icra_B->IntegratorLimited =
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE;

    /* Integrator: '<S3>/Integrator Limited1' */
    /* Limited  Integrator  */
    if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE >= max_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE = max_u;
    } else if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE <=
               min_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE = min_u;
    }

    /* Integrator: '<S3>/Integrator Limited1' */
    uavNavigationRL2D_2022a_icra_B->IntegratorLimited1 =
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE;

    /* Integrator: '<S3>/Integrator Limited2' */
    /* Limited  Integrator  */
    if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE >= max_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE = max_u;
    } else if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE <=
               min_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE = min_u;
    }

    /* Integrator: '<S3>/Integrator Limited2' */
    uavNavigationRL2D_2022a_icra_B->IntegratorLimited2 =
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE;

    /* Integrator: '<S3>/Integrator Limited3' */
    /* Limited  Integrator  */
    if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE >= max_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE = max_u;
    } else if (uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE <=
               min_u) {
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE = min_u;
    }

    /* Integrator: '<S3>/Integrator Limited3' */
    uavNavigationRL2D_2022a_icra_B->IntegratorLimited3 =
      uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE;

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[0];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[0];
    if (uavNavigationRL2D_2022a_icra_B->IntegratorLimited > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->IntegratorLimited < q_idx_1)) {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->IntegratorLimited;
    }

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[0] = q_idx_1;

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[1];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[1];
    if (uavNavigationRL2D_2022a_icra_B->IntegratorLimited1 > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->IntegratorLimited1 < q_idx_1))
    {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->IntegratorLimited1;
    }

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[1] = q_idx_1;

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[2];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[2];
    if (uavNavigationRL2D_2022a_icra_B->IntegratorLimited2 > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->IntegratorLimited2 < q_idx_1))
    {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->IntegratorLimited2;
    }

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[2] = q_idx_1;

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[3];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[3];
    if (uavNavigationRL2D_2022a_icra_B->IntegratorLimited3 > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->IntegratorLimited3 < q_idx_1))
    {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->IntegratorLimited3;
    }

    /* Saturate: '<S5>/Maximum  Command Authority1' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[3] = q_idx_1;

    /* Product: '<S5>/Product1' incorporates:
     *  Product: '<S5>/Matrix Divide'
     */
    tmp_2 = &uavNavigationRL2D_2022a_icra_B->MatrixDivide[0];
    b_y[0] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[0];
    b_y[1] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[1];
    b_y[2] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[2];
    b_y[3] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[3];
    for (i = 0; i <= 2; i += 2) {
      /* Product: '<S5>/Product1' */
      tmp_0 = _mm_loadu_pd(&tmp_2[i]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(b_y[0]));
      tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&tmp_2[i + 4]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(b_y[1]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S5>/Product1' */
      tmp_1 = _mm_loadu_pd(&tmp_2[i + 8]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(b_y[2]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S5>/Product1' */
      tmp_1 = _mm_loadu_pd(&tmp_2[i + 12]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(b_y[3]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S5>/Product1' */
      _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->Product1[i], tmp_0);
    }

    /* Gain: '<S5>/Gain13' */
    uavNavigationRL2D_2022a_icra_B->Gain13 =
      uavNavigationRL2D_2022a_icra_P.Gain13_Gain *
      uavNavigationRL2D_2022a_icra_B->Product1[3];

    /* Clock: '<S41>/Clock' */
    uavNavigationRL2D_2022a_icra_B->Clock =
      uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      boolean_T mask2;

      /* Memory: '<S41>/Memory' */
      uavNavigationRL2D_2022a_icra_B->Memory[0] =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[0];
      uavNavigationRL2D_2022a_icra_B->Memory[1] =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[1];
      uavNavigationRL2D_2022a_icra_B->Memory[2] =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[2];
      uavNavigationRL2D_2022a_icra_B->Memory[3] =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[3];

      /* MATLAB Function: '<S41>/kinematic differential equation' */
      tmp[0] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[1];
      tmp[4] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[2];
      tmp[8] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[3];
      tmp[1] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[0];
      tmp[5] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[3];
      tmp[9] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[2];
      tmp[2] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[3];
      tmp[6] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[0];
      tmp[10] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[1];
      tmp[3] = 0.5 * -uavNavigationRL2D_2022a_icra_B->Memory[2];
      tmp[7] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[1];
      tmp[11] = 0.5 * uavNavigationRL2D_2022a_icra_B->Memory[0];
      st_idx_1 = uavNavigationRL2D_2022a_icra_B->TransferFcn7;
      q_idx_2 = uavNavigationRL2D_2022a_icra_B->TransferFcn8;
      q_idx_3 = uavNavigationRL2D_2022a_icra_B->Gain13;
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->Clock -
        uavNavigationRL2D_2022a_icra_DW->last_t;
      st_idx_2 = 0.0;
      scale = 3.3121686421112381E-170;
      for (i = 0; i < 4; i++) {
        uavNavigationRL2D_2022a_icra_B->quat_rate[i] = 0.0;
        uavNavigationRL2D_2022a_icra_B->quat_rate[i] += tmp[i] * st_idx_1;
        uavNavigationRL2D_2022a_icra_B->quat_rate[i] += tmp[i + 4] * q_idx_2;
        uavNavigationRL2D_2022a_icra_B->quat_rate[i] += tmp[i + 8] * q_idx_3;
        uavNavigationRL2D_2022a_icra_B->quat[i] =
          uavNavigationRL2D_2022a_icra_B->quat_rate[i] * q_idx_1 +
          uavNavigationRL2D_2022a_icra_B->Memory[i];
        absxk = std::abs(uavNavigationRL2D_2022a_icra_B->quat[i]);
        if (absxk > scale) {
          t = scale / absxk;
          st_idx_2 = st_idx_2 * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          st_idx_2 += t * t;
        }
      }

      st_idx_2 = scale * std::sqrt(st_idx_2);
      uavNavigationRL2D_2022a_icra_B->quat[0] /= st_idx_2;
      b_y[0] = uavNavigationRL2D_2022a_icra_B->quat[0] *
        uavNavigationRL2D_2022a_icra_B->quat[0];
      uavNavigationRL2D_2022a_icra_B->quat[1] /= st_idx_2;
      b_y[1] = uavNavigationRL2D_2022a_icra_B->quat[1] *
        uavNavigationRL2D_2022a_icra_B->quat[1];
      uavNavigationRL2D_2022a_icra_B->quat[2] /= st_idx_2;
      b_y[2] = uavNavigationRL2D_2022a_icra_B->quat[2] *
        uavNavigationRL2D_2022a_icra_B->quat[2];
      uavNavigationRL2D_2022a_icra_B->quat[3] /= st_idx_2;
      b_y[3] = uavNavigationRL2D_2022a_icra_B->quat[3] *
        uavNavigationRL2D_2022a_icra_B->quat[3];
      roll_pitch_DCM[0] = ((uavNavigationRL2D_2022a_icra_B->quat[0] *
                            uavNavigationRL2D_2022a_icra_B->quat[0] +
                            uavNavigationRL2D_2022a_icra_B->quat[1] *
                            uavNavigationRL2D_2022a_icra_B->quat[1]) -
                           uavNavigationRL2D_2022a_icra_B->quat[2] *
                           uavNavigationRL2D_2022a_icra_B->quat[2]) -
        uavNavigationRL2D_2022a_icra_B->quat[3] *
        uavNavigationRL2D_2022a_icra_B->quat[3];
      roll_pitch_DCM[3] = (uavNavigationRL2D_2022a_icra_B->quat[1] *
                           uavNavigationRL2D_2022a_icra_B->quat[2] +
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[3]) * 2.0;
      roll_pitch_DCM[6] = (uavNavigationRL2D_2022a_icra_B->quat[1] *
                           uavNavigationRL2D_2022a_icra_B->quat[3] -
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[2]) * 2.0;
      roll_pitch_DCM[1] = (uavNavigationRL2D_2022a_icra_B->quat[1] *
                           uavNavigationRL2D_2022a_icra_B->quat[2] -
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[3]) * 2.0;
      roll_pitch_DCM[4] = ((uavNavigationRL2D_2022a_icra_B->quat[0] *
                            uavNavigationRL2D_2022a_icra_B->quat[0] -
                            uavNavigationRL2D_2022a_icra_B->quat[1] *
                            uavNavigationRL2D_2022a_icra_B->quat[1]) +
                           uavNavigationRL2D_2022a_icra_B->quat[2] *
                           uavNavigationRL2D_2022a_icra_B->quat[2]) -
        uavNavigationRL2D_2022a_icra_B->quat[3] *
        uavNavigationRL2D_2022a_icra_B->quat[3];
      roll_pitch_DCM[7] = (uavNavigationRL2D_2022a_icra_B->quat[2] *
                           uavNavigationRL2D_2022a_icra_B->quat[3] +
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[1]) * 2.0;
      roll_pitch_DCM[2] = (uavNavigationRL2D_2022a_icra_B->quat[1] *
                           uavNavigationRL2D_2022a_icra_B->quat[3] +
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[2]) * 2.0;
      roll_pitch_DCM[5] = (uavNavigationRL2D_2022a_icra_B->quat[2] *
                           uavNavigationRL2D_2022a_icra_B->quat[3] -
                           uavNavigationRL2D_2022a_icra_B->quat[0] *
                           uavNavigationRL2D_2022a_icra_B->quat[1]) * 2.0;
      roll_pitch_DCM[8] = ((uavNavigationRL2D_2022a_icra_B->quat[0] *
                            uavNavigationRL2D_2022a_icra_B->quat[0] -
                            uavNavigationRL2D_2022a_icra_B->quat[1] *
                            uavNavigationRL2D_2022a_icra_B->quat[1]) -
                           uavNavigationRL2D_2022a_icra_B->quat[2] *
                           uavNavigationRL2D_2022a_icra_B->quat[2]) +
        uavNavigationRL2D_2022a_icra_B->quat[3] *
        uavNavigationRL2D_2022a_icra_B->quat[3];
      st_idx_2 = b_y[0];
      st_idx_2 += b_y[1];
      st_idx_2 += b_y[2];
      st_idx_2 += b_y[3];
      st_idx_2 = 1.0 / std::sqrt(st_idx_2);
      b_y[0] = uavNavigationRL2D_2022a_icra_B->quat[0] * st_idx_2;
      b_y[1] = uavNavigationRL2D_2022a_icra_B->quat[1] * st_idx_2;
      b_y[2] = uavNavigationRL2D_2022a_icra_B->quat[2] * st_idx_2;
      b_y[3] = uavNavigationRL2D_2022a_icra_B->quat[3] * st_idx_2;
      st_idx_2 = (b_y[1] * b_y[3] - b_y[0] * b_y[2]) * -2.0;
      mask1 = (st_idx_2 >= 0.99999999999999778);
      mask2 = (st_idx_2 <= -0.99999999999999778);
      if (mask1) {
        st_idx_2 = 1.0;
      }

      if (mask2) {
        st_idx_2 = -1.0;
      }

      mask1 = (mask1 || mask2);
      q_idx_2 = rt_atan2d_snf((b_y[1] * b_y[2] + b_y[0] * b_y[3]) * 2.0, ((b_y[0]
        * b_y[0] + b_y[1] * b_y[1]) - b_y[2] * b_y[2]) - b_y[3] * b_y[3]);
      q_idx_3 = std::asin(st_idx_2);
      q_idx_1 = rt_atan2d_snf((b_y[2] * b_y[3] + b_y[0] * b_y[1]) * 2.0, ((b_y[0]
        * b_y[0] - b_y[1] * b_y[1]) - b_y[2] * b_y[2]) + b_y[3] * b_y[3]);
      euler_angle[0] = q_idx_2;
      euler_angle[1] = q_idx_3;
      euler_angle[2] = q_idx_1;
      i = 0;
      if (mask1) {
        i = 1;
      }

      x_size = i;
      if (i - 1 >= 0) {
        x_data = st_idx_2;
      }

      for (iy = 0; iy < i; iy++) {
        q_idx_1 = x_data;
        if (std::isnan(q_idx_1)) {
          q_idx_1 = (rtNaN);
        } else if (q_idx_1 < 0.0) {
          q_idx_1 = -1.0;
        } else {
          q_idx_1 = (q_idx_1 > 0.0);
        }

        x_data = q_idx_1;
      }

      i = 0;
      if (mask1) {
        i = 1;
      }

      f_size[0] = 1;
      f_size[1] = i;
      if (mask1) {
        g_data = 1;
      }

      i = f_size[1];
      if (x_size == i) {
        /* MATLAB Function: '<S41>/kinematic differential equation' */
        i = x_size;
        if (i - 1 >= 0) {
          q_idx_1 = b_y[1];
          st_idx_2 = b_y[0];
          q_idx_1 = rt_atan2d_snf(q_idx_1, st_idx_2);
          euler_angle[0] = -x_data * 2.0 * q_idx_1;
        }
      } else {
        /* MATLAB Function: '<S41>/kinematic differential equation' */
        uavNavigationR_binary_expand_op(euler_angle, &g_data, &x_data, &x_size,
          b_y, f_size);
      }

      /* MATLAB Function: '<S41>/kinematic differential equation' */
      st_idx_2 = euler_angle[0];
      q_idx_1 = st_idx_2;
      q_idx_3 = q_idx_1;
      q_idx_3 = std::cos(q_idx_3);
      q_idx_1 = std::sin(q_idx_1);
      st_idx_2 = q_idx_1;
      ct[2] = q_idx_3;
      a[0] = ct[2];
      a[3] = -st_idx_2;
      a[6] = 0.0;
      a[1] = 0.0 * ct[2] * 0.0 + st_idx_2;
      a[4] = ct[2] - 0.0 * st_idx_2;
      a[7] = -0.0;
      a[2] = 0.0 * st_idx_2 - ct[2] * 0.0;
      a[5] = 0.0 * ct[2] + 0.0 * st_idx_2;
      a[8] = 1.0;
      for (i = 0; i < 3; i++) {
        for (iy = 0; iy < 3; iy++) {
          tempR[i + 3 * iy] = 0.0;
          tempR[i + 3 * iy] += a[3 * iy] * roll_pitch_DCM[3 * i];
          tempR[i + 3 * iy] += a[3 * iy + 1] * roll_pitch_DCM[3 * i + 1];
          tempR[i + 3 * iy] += a[3 * iy + 2] * roll_pitch_DCM[3 * i + 2];
        }
      }

      mask1 = true;
      for (iy = 0; iy < 9; iy++) {
        st_idx_2 = tempR[iy];
        if (mask1 && ((!std::isinf(st_idx_2)) && (!std::isnan(st_idx_2)))) {
        } else {
          mask1 = false;
        }
      }

      if (mask1) {
        uavNavigationRL2D_2022a_icr_svd(tempR, d_y, ct, d_b);
      } else {
        for (i = 0; i < 9; i++) {
          d_y[i] = (rtNaN);
          d_b[i] = (rtNaN);
        }
      }

      for (i = 0; i < 3; i++) {
        for (iy = 0; iy < 3; iy++) {
          a[i + 3 * iy] = 0.0;
          q_idx_1 = a[3 * iy + i];
          q_idx_1 += d_y[i] * d_b[iy];
          a[i + 3 * iy] = q_idx_1;
          q_idx_1 = a[3 * iy + i];
          q_idx_1 += d_y[i + 3] * d_b[iy + 3];
          a[i + 3 * iy] = q_idx_1;
          q_idx_1 = a[3 * iy + i];
          q_idx_1 += d_y[i + 6] * d_b[iy + 6];
          a[i + 3 * iy] = q_idx_1;
        }
      }

      t = a[0];
      t += a[4];
      t += a[8];
      st_idx_2 = std::acos((t - 1.0) / 2.0);
      ct[0] = a[5] - a[7];
      ct[1] = a[6] - a[2];
      if (std::sin(st_idx_2) >= 0.0001) {
        q_idx_1 = 1.0 / (2.0 * std::sin(st_idx_2));
        q_idx_3 = ct[0];
        q_idx_3 = q_idx_3 * q_idx_1 * st_idx_2;
        ct[0] = q_idx_3;
        q_idx_3 = ct[1];
        q_idx_3 = q_idx_3 * q_idx_1 * st_idx_2;
        ct[1] = q_idx_3;
      } else if (t - 1.0 > 0.0) {
        st_idx_2 = 0.5 - (t - 3.0) / 12.0;
        q_idx_3 = ct[0];
        q_idx_3 *= st_idx_2;
        ct[0] = q_idx_3;
        q_idx_3 = ct[1];
        q_idx_3 *= st_idx_2;
        ct[1] = q_idx_3;
      } else {
        ct[0] = a[0];
        ct[1] = a[4];
        ct[2] = a[8];
        if (!std::isnan(ct[0])) {
          iy = 0;
        } else {
          iy = -1;
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i < 4)) {
            if (!std::isnan(ct[i - 1])) {
              iy = i - 1;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (iy + 1 == 0) {
          x_size = 0;
        } else {
          scale = ct[iy];
          x_size = iy;
          for (i = iy + 2; i < 4; i++) {
            if (scale < ct[i - 1]) {
              scale = ct[i - 1];
              x_size = i - 1;
            }
          }
        }

        scale = std::fmod(static_cast<real_T>(x_size) + 1.0, 3.0);
        absxk = std::fmod((static_cast<real_T>(x_size) + 1.0) + 1.0, 3.0);
        t = std::sqrt(((a[3 * x_size + x_size] - a[3 * static_cast<int32_T>
                        (scale) + static_cast<int32_T>(scale)]) - a[3 *
                       static_cast<int32_T>(absxk) + static_cast<int32_T>(absxk)])
                      + 1.0);
        ct[0] = 0.0;
        ct[1] = 0.0;
        ct[2] = 0.0;
        ct[x_size] = t / 2.0;
        ct[static_cast<int32_T>(scale)] = (a[3 * x_size + static_cast<int32_T>
          (scale)] + a[3 * static_cast<int32_T>(scale) + x_size]) / (2.0 * t);
        ct[static_cast<int32_T>(absxk)] = (a[3 * x_size + static_cast<int32_T>
          (absxk)] + a[3 * static_cast<int32_T>(absxk) + x_size]) / (2.0 * t);
        scale = 3.3121686421112381E-170;
        absxk = std::abs(ct[0]);
        if (absxk > 3.3121686421112381E-170) {
          q_idx_1 = 1.0;
          scale = absxk;
        } else {
          t = absxk / 3.3121686421112381E-170;
          q_idx_1 = t * t;
        }

        absxk = std::abs(ct[1]);
        if (absxk > scale) {
          t = scale / absxk;
          q_idx_1 = q_idx_1 * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          q_idx_1 += t * t;
        }

        absxk = std::abs(ct[2]);
        if (absxk > scale) {
          t = scale / absxk;
          q_idx_1 = q_idx_1 * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          q_idx_1 += t * t;
        }

        q_idx_1 = scale * std::sqrt(q_idx_1);
        q_idx_3 = ct[0];
        q_idx_3 = st_idx_2 * q_idx_3 / q_idx_1;
        ct[0] = q_idx_3;
        q_idx_3 = ct[1];
        q_idx_3 = st_idx_2 * q_idx_3 / q_idx_1;
        ct[1] = q_idx_3;
      }

      uavNavigationRL2D_2022a_icra_B->roll_h = ct[0];
      uavNavigationRL2D_2022a_icra_B->pitch = ct[1];
      for (i = 0; i < 3; i++) {
        uavNavigationRL2D_2022a_icra_B->DCM[3 * i] = roll_pitch_DCM[i];
        uavNavigationRL2D_2022a_icra_B->DCM[3 * i + 1] = roll_pitch_DCM[i + 3];
        uavNavigationRL2D_2022a_icra_B->DCM[3 * i + 2] = roll_pitch_DCM[i + 6];
      }

      uavNavigationRL2D_2022a_icra_DW->last_t =
        uavNavigationRL2D_2022a_icra_B->Clock;
      uavNavigationRL2D_2022a_icra_B->yaw = euler_angle[0];
    }

    /* TransportDelay: '<S5>/Transport Delay14' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay14_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime -
        uavNavigationRL2D_2022a_icra_P.TransportDelay14_Delay;
      if (uavNavigationRL2D_2022a_icra_P.TransportDelay14_Delay == 0.0)
        uavNavigationRL2D_2022a_icra_B->TransportDelay14 =
          uavNavigationRL2D_2022a_icra_B->yaw;
      else
        uavNavigationRL2D_2022a_icra_B->TransportDelay14 = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *uBuffer,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.CircularBufSize,
          &uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Last,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head,
          uavNavigationRL2D_2022a_icra_P.TransportDelay14_InitOutput,
          1,
          0);
    }

    /* Outport: '<Root>/states_output' */
    uavNavigationRL2D_states_output[0] = uavNavigationRL2D_2022a_icra_B->Sum1[0];
    uavNavigationRL2D_states_output[1] = uavNavigationRL2D_2022a_icra_B->Sum1[1];
    uavNavigationRL2D_states_output[2] = uavNavigationRL2D_2022a_icra_B->Sum1[2];
    uavNavigationRL2D_states_output[3] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay;
    uavNavigationRL2D_states_output[4] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay1;
    uavNavigationRL2D_states_output[5] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay2;
    uavNavigationRL2D_states_output[6] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay3;
    uavNavigationRL2D_states_output[7] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay4;
    uavNavigationRL2D_states_output[8] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay5;
    uavNavigationRL2D_states_output[9] = uavNavigationRL2D_2022a_icra_B->roll;
    uavNavigationRL2D_states_output[10] = uavNavigationRL2D_2022a_icra_B->Gain_n;
    uavNavigationRL2D_states_output[11] =
      uavNavigationRL2D_2022a_icra_B->TransportDelay14;

    /* MATLAB Function: '<Root>/MATLAB Function' incorporates:
     *  Constant: '<Root>/lower_limit_x'
     *  Constant: '<Root>/lower_limit_y'
     *  Constant: '<Root>/lower_limit_z'
     *  Constant: '<Root>/upper_limit_x'
     *  Constant: '<Root>/upper_limit_y'
     *  Constant: '<Root>/upper_limit_z'
     *  Inport: '<Root>/velocity_cmd_x'
     *  Inport: '<Root>/velocity_cmd_y'
     *  Inport: '<Root>/velocity_cmd_z'
     */
    uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim =
      uavNavigationRL2_velocity_cmd_x;
    uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim =
      uavNavigationRL2_velocity_cmd_y;
    uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim =
      uavNavigationRL2_velocity_cmd_z;
    if (uavNavigationRL2D_2022a_icra_B->Sum1[0] >
        uavNavigationRL2D_2022a_icra_P.upper_limit_x_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim = std::fmin
        (uavNavigationRL2_velocity_cmd_x, 0.0);
    }

    if (uavNavigationRL2D_2022a_icra_B->Sum1[0] <
        uavNavigationRL2D_2022a_icra_P.lower_limit_x_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim = std::fmax
        (uavNavigationRL2_velocity_cmd_x, 0.0);
    }

    if (uavNavigationRL2D_2022a_icra_B->Sum1[1] >
        uavNavigationRL2D_2022a_icra_P.upper_limit_y_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim = std::fmin
        (uavNavigationRL2_velocity_cmd_y, 0.0);
    }

    if (uavNavigationRL2D_2022a_icra_B->Sum1[1] <
        uavNavigationRL2D_2022a_icra_P.lower_limit_y_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim = std::fmax
        (uavNavigationRL2_velocity_cmd_y, 0.0);
    }

    if (uavNavigationRL2D_2022a_icra_B->Sum1[2] >
        uavNavigationRL2D_2022a_icra_P.upper_limit_z_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim = std::fmin
        (uavNavigationRL2_velocity_cmd_z, 0.0);
    }

    if (uavNavigationRL2D_2022a_icra_B->Sum1[2] <
        uavNavigationRL2D_2022a_icra_P.lower_limit_z_Value) {
      uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim = std::fmax
        (uavNavigationRL2_velocity_cmd_z, 0.0);
    }

    /* End of MATLAB Function: '<Root>/MATLAB Function' */

    /* Outport: '<Root>/velocity_references' */
    uavNavigati_velocity_references[0] =
      uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim;
    uavNavigati_velocity_references[1] =
      uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim;
    uavNavigati_velocity_references[2] =
      uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim;

    /* ZeroOrderHold: '<S4>/Zero-Order Hold3' incorporates:
     *  Constant: '<S2>/Constant1'
     *  MATLAB Function: '<S18>/reference correction'
     */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold3' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3 =
        uavNavigationRL2D_2022a_icra_B->TransportDelay14;

      /* Gain: '<S10>/Gain' */
      uavNavigationRL2D_2022a_icra_B->Gain_i =
        uavNavigationRL2D_2022a_icra_P.Gain_Gain_f *
        uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3;
      uavNavigationRL2D_2022a_icra_B->Constant1 =
        uavNavigationRL2D_2022a_icra_P.Constant1_Value_d;

      /* MATLAB Function: '<S18>/reference correction' incorporates:
       *  Constant: '<S2>/Constant1'
       */
      st_idx_1 = std::cos(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_3 = std::sin(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_2 = std::sin(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_1 = std::cos(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_2 = -q_idx_2;
      scale = uavNavigationRL2D_2022a_icra_B->Constant1;
      absxk = uavNavigationRL2D_2022a_icra_B->Constant1;
      st_idx_2 = st_idx_1 * scale;
      out[0] = st_idx_2;
      st_idx_2 = out[0];

      /* MATLAB Function: '<S18>/reference correction' */
      st_idx_2 += q_idx_3 * absxk;
      out[0] = st_idx_2;

      /* MATLAB Function: '<S18>/reference correction' */
      st_idx_2 = q_idx_2 * scale;
      out[1] = st_idx_2;
      st_idx_2 = out[1];

      /* MATLAB Function: '<S18>/reference correction' */
      st_idx_2 += q_idx_1 * absxk;
      out[1] = st_idx_2;

      /* MATLAB Function: '<S18>/reference correction' */
      uavNavigationRL2D_2022a_icra_B->x_out = out[0];
      uavNavigationRL2D_2022a_icra_B->y_out = out[1];

      /* ZeroOrderHold: '<S4>/Zero-Order Hold8' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold8 =
        uavNavigationRL2D_2022a_icra_B->Sum1[0];

      /* ZeroOrderHold: '<S4>/Zero-Order Hold9' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold9 =
        uavNavigationRL2D_2022a_icra_B->Sum1[1];

      /* MATLAB Function: '<S19>/reference correction' */
      uavNavigati_referencecorrection
        (uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3,
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold8,
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold9,
         &uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_c);

      /* Sum: '<S9>/Sum3' */
      uavNavigationRL2D_2022a_icra_B->Sum3 =
        uavNavigationRL2D_2022a_icra_B->x_out -
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_c.x_out;

      /* Product: '<S9>/Product' incorporates:
       *  Constant: '<Root>/True1'
       */
      uavNavigationRL2D_2022a_icra_B->Product =
        uavNavigationRL2D_2022a_icra_P.True1_Value *
        uavNavigationRL2D_2022a_icra_B->Sum3;

      /* Gain: '<S22>/Gain5' */
      uavNavigationRL2D_2022a_icra_B->Gain5 =
        uavNavigationRL2D_2022a_icra_P.Gain5_Gain *
        uavNavigationRL2D_2022a_icra_B->Product;
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold3' */

    /* TransferFcn: '<S22>/Transfer Fcn1' */
    uavNavigationRL2D_2022a_icra_B->TransferFcn1 =
      uavNavigationRL2D_2022a_icra_P.TransferFcn1_C *
      uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE;
    uavNavigationRL2D_2022a_icra_B->TransferFcn1 +=
      uavNavigationRL2D_2022a_icra_P.TransferFcn1_D *
      uavNavigationRL2D_2022a_icra_B->Gain5;

    /* RelationalOperator: '<S22>/Relational Operator' incorporates:
     *  Constant: '<S22>/Constant4'
     */
    uavNavigationRL2D_2022a_icra_B->RelationalOperator =
      (uavNavigationRL2D_2022a_icra_B->TransferFcn1 >
       uavNavigationRL2D_2022a_icra_P.Constant4_Value_a);
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S22>/Memory2' */
      uavNavigationRL2D_2022a_icra_B->Memory2 =
        uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput;

      /* Memory: '<S23>/Memory' */
      uavNavigationRL2D_2022a_icra_B->Memory_k =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_l;
    }

    /* Logic: '<S22>/Logical Operator1' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator1 =
      uavNavigationRL2D_2022a_icra_B->RelationalOperator ^
      uavNavigationRL2D_2022a_icra_B->Memory2;

    /* Logic: '<S22>/Logical Operator' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator =
      !uavNavigationRL2D_2022a_icra_B->RelationalOperator;

    /* Logic: '<S22>/Logical Operator4' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator4 =
      (uavNavigationRL2D_2022a_icra_B->LogicalOperator1 &&
       uavNavigationRL2D_2022a_icra_B->LogicalOperator);

    /* Switch: '<S23>/Switch' */
    if (uavNavigationRL2D_2022a_icra_B->LogicalOperator4) {
      /* Switch: '<S23>/Switch' incorporates:
       *  Constant: '<S23>/Constant2'
       */
      uavNavigationRL2D_2022a_icra_B->Switch =
        uavNavigationRL2D_2022a_icra_P.Constant2_Value_e;
    } else {
      /* Switch: '<S23>/Switch' */
      uavNavigationRL2D_2022a_icra_B->Switch =
        uavNavigationRL2D_2022a_icra_B->Memory_k;
    }

    /* End of Switch: '<S23>/Switch' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* Gain: '<S22>/Gain4' */
      uavNavigationRL2D_2022a_icra_B->Gain4 =
        uavNavigationRL2D_2022a_icra_P.Gain4_Gain *
        uavNavigationRL2D_2022a_icra_B->Product;
    }

    /* MinMax: '<S23>/MinMax' */
    st_idx_2 = std::fmax(uavNavigationRL2D_2022a_icra_B->Switch,
                         uavNavigationRL2D_2022a_icra_B->Gain4);

    /* MinMax: '<S23>/MinMax' */
    uavNavigationRL2D_2022a_icra_B->MinMax = st_idx_2;

    /* Logic: '<S22>/Logical Operator3' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator3 =
      (uavNavigationRL2D_2022a_icra_B->LogicalOperator1 &&
       uavNavigationRL2D_2022a_icra_B->RelationalOperator);
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S24>/Memory1' */
      uavNavigationRL2D_2022a_icra_B->Memory1 =
        uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput;
    }

    /* Switch: '<S24>/Switch1' */
    if (uavNavigationRL2D_2022a_icra_B->LogicalOperator3) {
      /* Switch: '<S24>/Switch1' incorporates:
       *  Constant: '<S24>/Constant3'
       */
      uavNavigationRL2D_2022a_icra_B->Switch1 =
        uavNavigationRL2D_2022a_icra_P.Constant3_Value;
    } else {
      /* Switch: '<S24>/Switch1' */
      uavNavigationRL2D_2022a_icra_B->Switch1 =
        uavNavigationRL2D_2022a_icra_B->Memory1;
    }

    /* End of Switch: '<S24>/Switch1' */

    /* MinMax: '<S24>/MinMax1' */
    st_idx_2 = std::fmin(uavNavigationRL2D_2022a_icra_B->Gain4,
                         uavNavigationRL2D_2022a_icra_B->Switch1);

    /* MinMax: '<S24>/MinMax1' */
    uavNavigationRL2D_2022a_icra_B->MinMax1 = st_idx_2;
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Switch: '<S15>/Switch1' */
      if (uavNavigationRL2D_2022a_icra_P.Switch1_Threshold < 0.0) {
        /* Switch: '<S15>/Switch1' incorporates:
         *  Constant: '<S15>/Constant'
         */
        uavNavigationRL2D_2022a_icra_B->Switch1_e =
          uavNavigationRL2D_2022a_icra_P.Constant_Value;
      } else {
        /* Switch: '<S15>/Switch1' incorporates:
         *  Constant: '<S15>/Constant1'
         */
        uavNavigationRL2D_2022a_icra_B->Switch1_e =
          uavNavigationRL2D_2022a_icra_P.Constant1_Value;
      }

      /* End of Switch: '<S15>/Switch1' */
    }

    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* Product: '<S15>/Product2' incorporates:
       *  Constant: '<S15>/Constant2'
       */
      uavNavigationRL2D_2022a_icra_B->Product2 = kp_x *
        uavNavigationRL2D_2022a_icra_B->Product;

      /* ZeroOrderHold: '<S4>/Zero-Order Hold2' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[0] =
        uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[1] =
        uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[2] =
        uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim;

      /* MATLAB Function: '<S21>/reference correction' */
      uavNavigati_referencecorrection
        (uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3,
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[0],
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[1],
         &uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_l);
    }

    /* Derivative: '<S5>/Derivative3' */
    if ((uavNavigationRL2D_2022a_icra_DW->TimeStampA >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0]) &&
        (uavNavigationRL2D_2022a_icra_DW->TimeStampB >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0])) {
      /* Derivative: '<S5>/Derivative3' */
      uavNavigationRL2D_2022a_icra_B->Derivative3 = 0.0;
    } else {
      st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampA;
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA;
      if (uavNavigationRL2D_2022a_icra_DW->TimeStampA <
          uavNavigationRL2D_2022a_icra_DW->TimeStampB) {
        if (uavNavigationRL2D_2022a_icra_DW->TimeStampB <
            uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
          st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB;
          lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB;
        }
      } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA >=
                 uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
        st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB;
        lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB;
      }

      st_idx_2 = uavNavigationRL2D_2022a_icra_M->Timing.t[0] - st_idx_2;

      /* Derivative: '<S5>/Derivative3' */
      uavNavigationRL2D_2022a_icra_B->Derivative3 =
        (uavNavigationRL2D_2022a_icra_B->Sum1[0] - *lastU) / st_idx_2;
    }

    /* End of Derivative: '<S5>/Derivative3' */

    /* ZeroOrderHold: '<S4>/Zero-Order Hold10' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold10' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold10 =
        uavNavigationRL2D_2022a_icra_B->Derivative3;
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold10' */

    /* Derivative: '<S5>/Derivative4' */
    if ((uavNavigationRL2D_2022a_icra_DW->TimeStampA_b >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0]) &&
        (uavNavigationRL2D_2022a_icra_DW->TimeStampB_l >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0])) {
      /* Derivative: '<S5>/Derivative4' */
      uavNavigationRL2D_2022a_icra_B->Derivative4 = 0.0;
    } else {
      st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampA_b;
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_i;
      if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_b <
          uavNavigationRL2D_2022a_icra_DW->TimeStampB_l) {
        if (uavNavigationRL2D_2022a_icra_DW->TimeStampB_l <
            uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
          st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB_l;
          lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_h;
        }
      } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_b >=
                 uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
        st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB_l;
        lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_h;
      }

      st_idx_2 = uavNavigationRL2D_2022a_icra_M->Timing.t[0] - st_idx_2;

      /* Derivative: '<S5>/Derivative4' */
      uavNavigationRL2D_2022a_icra_B->Derivative4 =
        (uavNavigationRL2D_2022a_icra_B->Sum1[1] - *lastU) / st_idx_2;
    }

    /* End of Derivative: '<S5>/Derivative4' */

    /* ZeroOrderHold: '<S4>/Zero-Order Hold11' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold11' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold11 =
        uavNavigationRL2D_2022a_icra_B->Derivative4;

      /* MATLAB Function: '<S20>/reference correction' */
      uavNavigati_referencecorrection
        (uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3,
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold10,
         uavNavigationRL2D_2022a_icra_B->ZeroOrderHold11,
         &uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_n);

      /* Sum: '<S9>/Sum6' */
      uavNavigationRL2D_2022a_icra_B->Sum6 =
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_l.x_out -
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_n.x_out;

      /* Product: '<S15>/Product3' incorporates:
       *  Constant: '<S15>/Constant5'
       */
      uavNavigationRL2D_2022a_icra_B->Product3 = kd_x *
        uavNavigationRL2D_2022a_icra_B->Sum6;

      /* Sum: '<S15>/Sum1' */
      uavNavigationRL2D_2022a_icra_B->Sum1_f =
        uavNavigationRL2D_2022a_icra_B->Product2 +
        uavNavigationRL2D_2022a_icra_B->Product3;

      /* Sum: '<S9>/Sum5' */
      uavNavigationRL2D_2022a_icra_B->Sum5 =
        uavNavigationRL2D_2022a_icra_B->y_out -
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_c.y_out;

      /* Product: '<S9>/Product1' incorporates:
       *  Constant: '<Root>/True2'
       */
      uavNavigationRL2D_2022a_icra_B->Product1_o =
        uavNavigationRL2D_2022a_icra_P.True2_Value *
        uavNavigationRL2D_2022a_icra_B->Sum5;

      /* Product: '<S16>/Product2' incorporates:
       *  Constant: '<S16>/Constant3'
       */
      uavNavigationRL2D_2022a_icra_B->Product2_p = kp_y *
        uavNavigationRL2D_2022a_icra_B->Product1_o;

      /* Sum: '<S9>/Sum4' */
      uavNavigationRL2D_2022a_icra_B->Sum4 =
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_l.y_out -
        uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_n.y_out;

      /* Product: '<S16>/Product3' incorporates:
       *  Constant: '<S16>/Constant4'
       */
      uavNavigationRL2D_2022a_icra_B->Product3_o = kd_y *
        uavNavigationRL2D_2022a_icra_B->Sum4;

      /* Sum: '<S16>/Sum1' */
      uavNavigationRL2D_2022a_icra_B->Sum1_k =
        uavNavigationRL2D_2022a_icra_B->Product2_p +
        uavNavigationRL2D_2022a_icra_B->Product3_o;
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold11' */

    /* Switch: '<S15>/Switch' */
    if (uavNavigationRL2D_2022a_icra_B->Switch1_e >
        uavNavigationRL2D_2022a_icra_P.Switch_Threshold_f) {
      /* Switch: '<S22>/Switch2' */
      if (uavNavigationRL2D_2022a_icra_B->RelationalOperator) {
        /* Switch: '<S22>/Switch2' */
        uavNavigationRL2D_2022a_icra_B->Switch2 =
          uavNavigationRL2D_2022a_icra_B->MinMax;
      } else {
        /* Switch: '<S22>/Switch2' */
        uavNavigationRL2D_2022a_icra_B->Switch2 =
          uavNavigationRL2D_2022a_icra_B->MinMax1;
      }

      /* End of Switch: '<S22>/Switch2' */

      /* Gain: '<S22>/Beta' */
      uavNavigationRL2D_2022a_icra_B->Beta =
        uavNavigationRL2D_2022a_icra_P.Beta_Gain_m *
        uavNavigationRL2D_2022a_icra_B->Switch2;

      /* Sum: '<S22>/Sum4' */
      uavNavigationRL2D_2022a_icra_B->Sum4_l =
        uavNavigationRL2D_2022a_icra_B->Product +
        uavNavigationRL2D_2022a_icra_B->Beta;

      /* Signum: '<S22>/Sign' */
      st_idx_2 = uavNavigationRL2D_2022a_icra_B->Sum4_l;
      if (std::isnan(st_idx_2)) {
        /* Signum: '<S22>/Sign' */
        uavNavigationRL2D_2022a_icra_B->Sign = (rtNaN);
      } else if (st_idx_2 < 0.0) {
        /* Signum: '<S22>/Sign' */
        uavNavigationRL2D_2022a_icra_B->Sign = -1.0;
      } else {
        /* Signum: '<S22>/Sign' */
        uavNavigationRL2D_2022a_icra_B->Sign = (st_idx_2 > 0.0);
      }

      /* End of Signum: '<S22>/Sign' */

      /* Gain: '<S22>/h' */
      uavNavigationRL2D_2022a_icra_B->h =
        uavNavigationRL2D_2022a_icra_P.h_Gain_k *
        uavNavigationRL2D_2022a_icra_B->Sign;

      /* Switch: '<S15>/Switch' */
      uavNavigationRL2D_2022a_icra_B->Switch_b =
        uavNavigationRL2D_2022a_icra_B->h;
    } else {
      /* Switch: '<S15>/Switch' */
      uavNavigationRL2D_2022a_icra_B->Switch_b =
        uavNavigationRL2D_2022a_icra_B->Sum1_f;
    }

    /* End of Switch: '<S15>/Switch' */

    /* MATLAB Function: '<S34>/reference correction' */
    uavNavigati_referencecorrection(uavNavigationRL2D_2022a_icra_B->Gain_i,
      uavNavigationRL2D_2022a_icra_B->Switch_b,
      uavNavigationRL2D_2022a_icra_B->Sum1_k,
      &uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld);

    /* ZeroOrderHold: '<S4>/Zero-Order Hold6' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold6' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold6 =
        uavNavigationRL2D_2022a_icra_B->Sum1[2];

      /* Sum: '<S9>/Sum9' */
      uavNavigationRL2D_2022a_icra_B->Error =
        uavNavigationRL2D_2022a_icra_B->Constant1 -
        uavNavigationRL2D_2022a_icra_B->ZeroOrderHold6;

      /* Product: '<S9>/Product2' incorporates:
       *  Constant: '<Root>/False1'
       */
      uavNavigationRL2D_2022a_icra_B->Product2_a =
        uavNavigationRL2D_2022a_icra_P.False1_Value *
        uavNavigationRL2D_2022a_icra_B->Error;

      /* Gain: '<S25>/Gain5' */
      uavNavigationRL2D_2022a_icra_B->Gain5_h =
        uavNavigationRL2D_2022a_icra_P.Gain5_Gain_e *
        uavNavigationRL2D_2022a_icra_B->Product2_a;
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold6' */

    /* TransferFcn: '<S25>/Transfer Fcn1' */
    uavNavigationRL2D_2022a_icra_B->TransferFcn1_a =
      uavNavigationRL2D_2022a_icra_P.TransferFcn1_C_b *
      uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE_p;
    uavNavigationRL2D_2022a_icra_B->TransferFcn1_a +=
      uavNavigationRL2D_2022a_icra_P.TransferFcn1_D_d *
      uavNavigationRL2D_2022a_icra_B->Gain5_h;

    /* RelationalOperator: '<S25>/Relational Operator' incorporates:
     *  Constant: '<S25>/Constant4'
     */
    uavNavigationRL2D_2022a_icra_B->RelationalOperator_b =
      (uavNavigationRL2D_2022a_icra_B->TransferFcn1_a >
       uavNavigationRL2D_2022a_icra_P.Constant4_Value_f);
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S25>/Memory2' */
      uavNavigationRL2D_2022a_icra_B->Memory2_o =
        uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput_j;

      /* Memory: '<S26>/Memory' */
      uavNavigationRL2D_2022a_icra_B->Memory_g =
        uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_k;
    }

    /* Logic: '<S25>/Logical Operator1' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator1_m =
      uavNavigationRL2D_2022a_icra_B->RelationalOperator_b ^
      uavNavigationRL2D_2022a_icra_B->Memory2_o;

    /* Logic: '<S25>/Logical Operator' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator_c =
      !uavNavigationRL2D_2022a_icra_B->RelationalOperator_b;

    /* Logic: '<S25>/Logical Operator4' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator4_l =
      (uavNavigationRL2D_2022a_icra_B->LogicalOperator1_m &&
       uavNavigationRL2D_2022a_icra_B->LogicalOperator_c);

    /* Switch: '<S26>/Switch' */
    if (uavNavigationRL2D_2022a_icra_B->LogicalOperator4_l) {
      /* Switch: '<S26>/Switch' incorporates:
       *  Constant: '<S26>/Constant2'
       */
      uavNavigationRL2D_2022a_icra_B->Switch_l =
        uavNavigationRL2D_2022a_icra_P.Constant2_Value_f;
    } else {
      /* Switch: '<S26>/Switch' */
      uavNavigationRL2D_2022a_icra_B->Switch_l =
        uavNavigationRL2D_2022a_icra_B->Memory_g;
    }

    /* End of Switch: '<S26>/Switch' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* Gain: '<S25>/Gain4' */
      uavNavigationRL2D_2022a_icra_B->Gain4_c =
        uavNavigationRL2D_2022a_icra_P.Gain4_Gain_c *
        uavNavigationRL2D_2022a_icra_B->Product2_a;
    }

    /* MinMax: '<S26>/MinMax' */
    st_idx_2 = std::fmax(uavNavigationRL2D_2022a_icra_B->Switch_l,
                         uavNavigationRL2D_2022a_icra_B->Gain4_c);

    /* MinMax: '<S26>/MinMax' */
    uavNavigationRL2D_2022a_icra_B->MinMax_n = st_idx_2;

    /* Logic: '<S25>/Logical Operator3' */
    uavNavigationRL2D_2022a_icra_B->LogicalOperator3_f =
      (uavNavigationRL2D_2022a_icra_B->LogicalOperator1_m &&
       uavNavigationRL2D_2022a_icra_B->RelationalOperator_b);
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S27>/Memory1' */
      uavNavigationRL2D_2022a_icra_B->Memory1_e =
        uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput_o;
    }

    /* Switch: '<S27>/Switch1' */
    if (uavNavigationRL2D_2022a_icra_B->LogicalOperator3_f) {
      /* Switch: '<S27>/Switch1' incorporates:
       *  Constant: '<S27>/Constant3'
       */
      uavNavigationRL2D_2022a_icra_B->Switch1_j =
        uavNavigationRL2D_2022a_icra_P.Constant3_Value_a;
    } else {
      /* Switch: '<S27>/Switch1' */
      uavNavigationRL2D_2022a_icra_B->Switch1_j =
        uavNavigationRL2D_2022a_icra_B->Memory1_e;
    }

    /* End of Switch: '<S27>/Switch1' */

    /* MinMax: '<S27>/MinMax1' */
    st_idx_2 = std::fmin(uavNavigationRL2D_2022a_icra_B->Gain4_c,
                         uavNavigationRL2D_2022a_icra_B->Switch1_j);

    /* MinMax: '<S27>/MinMax1' */
    uavNavigationRL2D_2022a_icra_B->MinMax1_h = st_idx_2;
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Switch: '<S17>/Switch1' */
      if (uavNavigationRL2D_2022a_icra_P.Switch1_Threshold_n < 0.0) {
        /* Switch: '<S17>/Switch1' incorporates:
         *  Constant: '<S17>/Constant'
         */
        uavNavigationRL2D_2022a_icra_B->Switch1_g =
          uavNavigationRL2D_2022a_icra_P.Constant_Value_m;
      } else {
        /* Switch: '<S17>/Switch1' incorporates:
         *  Constant: '<S17>/Constant1'
         */
        uavNavigationRL2D_2022a_icra_B->Switch1_g =
          uavNavigationRL2D_2022a_icra_P.Constant1_Value_j;
      }

      /* End of Switch: '<S17>/Switch1' */
    }

    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* Product: '<S17>/Product2' incorporates:
       *  Constant: '<S17>/Constant3'
       */
      uavNavigationRL2D_2022a_icra_B->Product2_n = kp_z *
        uavNavigationRL2D_2022a_icra_B->Product2_a;
    }

    /* Derivative: '<S5>/Derivative5' */
    if ((uavNavigationRL2D_2022a_icra_DW->TimeStampA_h >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0]) &&
        (uavNavigationRL2D_2022a_icra_DW->TimeStampB_a >=
         uavNavigationRL2D_2022a_icra_M->Timing.t[0])) {
      /* Derivative: '<S5>/Derivative5' */
      uavNavigationRL2D_2022a_icra_B->Derivative5 = 0.0;
    } else {
      st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampA_h;
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_p;
      if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_h <
          uavNavigationRL2D_2022a_icra_DW->TimeStampB_a) {
        if (uavNavigationRL2D_2022a_icra_DW->TimeStampB_a <
            uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
          st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB_a;
          lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_o;
        }
      } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_h >=
                 uavNavigationRL2D_2022a_icra_M->Timing.t[0]) {
        st_idx_2 = uavNavigationRL2D_2022a_icra_DW->TimeStampB_a;
        lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_o;
      }

      st_idx_2 = uavNavigationRL2D_2022a_icra_M->Timing.t[0] - st_idx_2;

      /* Derivative: '<S5>/Derivative5' */
      uavNavigationRL2D_2022a_icra_B->Derivative5 =
        (uavNavigationRL2D_2022a_icra_B->Sum1[2] - *lastU) / st_idx_2;
    }

    /* End of Derivative: '<S5>/Derivative5' */

    /* ZeroOrderHold: '<S4>/Zero-Order Hold7' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold7' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold7 =
        uavNavigationRL2D_2022a_icra_B->Derivative5;

      /* Sum: '<S9>/Sum2' */
      uavNavigationRL2D_2022a_icra_B->Error_i =
        uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[2] -
        uavNavigationRL2D_2022a_icra_B->ZeroOrderHold7;

      /* Product: '<S17>/Product3' incorporates:
       *  Constant: '<S17>/Constant4'
       */
      uavNavigationRL2D_2022a_icra_B->Product3_p = kd_z *
        uavNavigationRL2D_2022a_icra_B->Error_i;

      /* Sum: '<S17>/Sum1' */
      uavNavigationRL2D_2022a_icra_B->Sum1_h =
        uavNavigationRL2D_2022a_icra_B->Product2_n +
        uavNavigationRL2D_2022a_icra_B->Product3_p;
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold7' */

    /* Switch: '<S9>/Switch' incorporates:
     *  Constant: '<Root>/False2'
     */
    if (uavNavigationRL2D_2022a_icra_P.False2_Value >
        uavNavigationRL2D_2022a_icra_P.Switch_Threshold_m) {
      /* Switch: '<S17>/Switch' */
      if (uavNavigationRL2D_2022a_icra_B->Switch1_g >
          uavNavigationRL2D_2022a_icra_P.Switch_Threshold) {
        /* Switch: '<S25>/Switch2' */
        if (uavNavigationRL2D_2022a_icra_B->RelationalOperator_b) {
          /* Switch: '<S25>/Switch2' */
          uavNavigationRL2D_2022a_icra_B->Switch2_f =
            uavNavigationRL2D_2022a_icra_B->MinMax_n;
        } else {
          /* Switch: '<S25>/Switch2' */
          uavNavigationRL2D_2022a_icra_B->Switch2_f =
            uavNavigationRL2D_2022a_icra_B->MinMax1_h;
        }

        /* End of Switch: '<S25>/Switch2' */

        /* Gain: '<S25>/Beta' */
        uavNavigationRL2D_2022a_icra_B->Beta_i =
          uavNavigationRL2D_2022a_icra_P.Beta_Gain *
          uavNavigationRL2D_2022a_icra_B->Switch2_f;

        /* Sum: '<S25>/Sum4' */
        uavNavigationRL2D_2022a_icra_B->Sum4_p =
          uavNavigationRL2D_2022a_icra_B->Product2_a +
          uavNavigationRL2D_2022a_icra_B->Beta_i;

        /* Signum: '<S25>/Sign' */
        st_idx_2 = uavNavigationRL2D_2022a_icra_B->Sum4_p;
        if (std::isnan(st_idx_2)) {
          /* Signum: '<S25>/Sign' */
          uavNavigationRL2D_2022a_icra_B->Sign_b = (rtNaN);
        } else if (st_idx_2 < 0.0) {
          /* Signum: '<S25>/Sign' */
          uavNavigationRL2D_2022a_icra_B->Sign_b = -1.0;
        } else {
          /* Signum: '<S25>/Sign' */
          uavNavigationRL2D_2022a_icra_B->Sign_b = (st_idx_2 > 0.0);
        }

        /* End of Signum: '<S25>/Sign' */

        /* Gain: '<S25>/h' */
        uavNavigationRL2D_2022a_icra_B->h_d =
          uavNavigationRL2D_2022a_icra_P.h_Gain *
          uavNavigationRL2D_2022a_icra_B->Sign_b;

        /* Switch: '<S17>/Switch' */
        uavNavigationRL2D_2022a_icra_B->Switch_p =
          uavNavigationRL2D_2022a_icra_B->h_d;
      } else {
        /* Switch: '<S17>/Switch' */
        uavNavigationRL2D_2022a_icra_B->Switch_p =
          uavNavigationRL2D_2022a_icra_B->Sum1_h;
      }

      /* End of Switch: '<S17>/Switch' */

      /* Switch: '<S9>/Switch' */
      uavNavigationRL2D_2022a_icra_B->Switch_h =
        uavNavigationRL2D_2022a_icra_B->Switch_p;
    } else {
      /* Switch: '<S9>/Switch' incorporates:
       *  Inport: '<Root>/u_cmd_z'
       */
      uavNavigationRL2D_2022a_icra_B->Switch_h = uavNavigationRL2D_2022a_u_cmd_z;
    }

    /* End of Switch: '<S9>/Switch' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S9>/Product3' incorporates:
       *  Constant: '<S5>/Constant'
       */
      st_idx_2 = 9.81 * T_aero / K_z;

      /* Product: '<S9>/Product3' incorporates:
       *  Constant: '<Root>/True3'
       */
      uavNavigationRL2D_2022a_icra_B->Product3_ox =
        uavNavigationRL2D_2022a_icra_P.True3_Value * st_idx_2;
    }

    /* Sum: '<S9>/Sum1' */
    uavNavigationRL2D_2022a_icra_B->Sum1_hp =
      uavNavigationRL2D_2022a_icra_B->Switch_h +
      uavNavigationRL2D_2022a_icra_B->Product3_ox;

    /* MATLAB Function: '<S33>/attitude mapper ' incorporates:
     *  Constant: '<Root>/True3'
     *  Constant: '<S2>/max_tilt_ref_limit_deg'
     *  Constant: '<S33>/virtual anti-gravity (just to make attitude defined)'
     */
    st_idx_2 = uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.x_out;
    q_idx_1 = uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.y_out;
    if ((uavNavigationRL2D_2022a_icra_P.True3_Value != 0.0) &&
        (uavNavigationRL2D_2022a_icra_P.max_tilt_ref_limit_deg_Value > 0.0)) {
      scale = std::tan(0.017453292519943295 *
                       uavNavigationRL2D_2022a_icra_P.max_tilt_ref_limit_deg_Value)
        * uavNavigationRL2D_2022a_icra_B->Sum1_hp;
      out[0] = uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.x_out;
      out[1] = uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.y_out;
      absxk = uavNavigationRL2D_2022a_ic_norm(out);
      if (absxk > scale) {
        st_idx_2 = out[0];
        st_idx_2 = st_idx_2 / absxk * scale;
        out[0] = st_idx_2;
        st_idx_2 = out[1];
        st_idx_2 = st_idx_2 / absxk * scale;
        out[1] = st_idx_2;
        st_idx_2 = out[0];
        q_idx_1 = out[1];
      }
    }

    euler_angle[0] = st_idx_2;
    euler_angle[1] = q_idx_1;
    euler_angle[2] = uavNavigationRL2D_2022a_icra_B->Sum1_hp;
    q_idx_1 = uavNavigationRL2D_2022a__norm_h(euler_angle);
    scale = euler_angle[0] / q_idx_1;
    st_idx_1 = euler_angle[1] / q_idx_1;
    st_idx_2 = euler_angle[2] / q_idx_1;
    if (!(uavNavigationRL2D_2022a_icra_P.True3_Value != 0.0)) {
      st_idx_2 = uavNavigationRL2D_2022a_icra_B->Sum1_hp;
      if (std::isnan(st_idx_2)) {
        st_idx_2 = (rtNaN);
      } else if (st_idx_2 < 0.0) {
        st_idx_2 = -1.0;
      } else {
        st_idx_2 = (st_idx_2 > 0.0);
      }

      q_idx_1 *= st_idx_2;
      euler_angle[0] = 0.0;
      euler_angle[1] = 0.0;
      euler_angle[2] =
        uavNavigationRL2D_2022a_icra_P.virtualantigravityjusttomakeatt;
      q_idx_2 = uavNavigationRL2D_2022a__norm_h(euler_angle);
      scale = 0.0 / q_idx_2;
      st_idx_1 = 0.0 / q_idx_2;
      st_idx_2 = euler_angle[2] / q_idx_2;
    }

    q_idx_2 = std::cos(uavNavigationRL2D_2022a_icra_B->Constant1);
    q_idx_3 = std::sin(uavNavigationRL2D_2022a_icra_B->Constant1);
    euler_angle[0] = q_idx_2;
    euler_angle[1] = q_idx_3;
    ct[0] = st_idx_1 * 0.0 - euler_angle[1] * st_idx_2;
    ct[1] = euler_angle[0] * st_idx_2 - scale * 0.0;
    ct[2] = scale * euler_angle[1] - euler_angle[0] * st_idx_1;
    q_idx_2 = uavNavigationRL2D_2022a__norm_h(ct);
    q_idx_3 = ct[0];
    q_idx_3 /= q_idx_2;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[3] = q_idx_3;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[6] = scale;
    ct[0] = q_idx_3;
    q_idx_3 = ct[1];
    q_idx_3 /= q_idx_2;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[4] = q_idx_3;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[7] = st_idx_1;
    ct[1] = q_idx_3;
    q_idx_3 = ct[2];
    q_idx_3 /= q_idx_2;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[5] = q_idx_3;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[8] = st_idx_2;
    ct[2] = q_idx_3;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[0] = ct[1] * st_idx_2 - st_idx_1 *
      ct[2];
    uavNavigationRL2D_2022a_icra_B->DCM_ref[1] = scale * ct[2] - ct[0] *
      st_idx_2;
    uavNavigationRL2D_2022a_icra_B->DCM_ref[2] = ct[0] * st_idx_1 - scale * ct[1];
    for (i = 0; i < 3; i++) {
      tempR[3 * i] = uavNavigationRL2D_2022a_icra_B->DCM_ref[i];
      tempR[3 * i + 1] = uavNavigationRL2D_2022a_icra_B->DCM_ref[i + 3];
      tempR[3 * i + 2] = uavNavigationRL2D_2022a_icra_B->DCM_ref[i + 6];
    }

    mask1 = true;
    for (iy = 0; iy < 9; iy++) {
      st_idx_2 = tempR[iy];
      if (mask1 && ((!std::isinf(st_idx_2)) && (!std::isnan(st_idx_2)))) {
      } else {
        mask1 = false;
      }
    }

    if (mask1) {
      uavNavigationRL2D_2022a_icr_svd(tempR, roll_pitch_DCM, euler_angle, a);
    } else {
      for (i = 0; i < 9; i++) {
        roll_pitch_DCM[i] = (rtNaN);
        a[i] = (rtNaN);
      }
    }

    for (i = 0; i < 3; i++) {
      for (iy = 0; iy < 3; iy++) {
        tempR[i + 3 * iy] = 0.0;
        tempR[i + 3 * iy] += roll_pitch_DCM[i] * a[iy];
        tempR[i + 3 * iy] += roll_pitch_DCM[i + 3] * a[iy + 3];
        tempR[i + 3 * iy] += roll_pitch_DCM[i + 6] * a[iy + 6];
      }
    }

    t = tempR[0];
    t += tempR[4];
    t += tempR[8];
    st_idx_2 = std::acos((t - 1.0) / 2.0);
    euler_angle[0] = tempR[5] - tempR[7];
    euler_angle[1] = tempR[6] - tempR[2];
    if (std::sin(st_idx_2) >= 0.0001) {
      scale = 1.0 / (2.0 * std::sin(st_idx_2));
      st_idx_1 = euler_angle[0];
      st_idx_1 = st_idx_1 * scale * st_idx_2;
      euler_angle[0] = st_idx_1;
      st_idx_1 = euler_angle[1];
      st_idx_1 = st_idx_1 * scale * st_idx_2;
      euler_angle[1] = st_idx_1;
    } else if (t - 1.0 > 0.0) {
      st_idx_2 = 0.5 - (t - 3.0) / 12.0;
      st_idx_1 = euler_angle[0];
      st_idx_1 *= st_idx_2;
      euler_angle[0] = st_idx_1;
      st_idx_1 = euler_angle[1];
      st_idx_1 *= st_idx_2;
      euler_angle[1] = st_idx_1;
    } else {
      euler_angle[0] = tempR[0];
      euler_angle[1] = tempR[4];
      euler_angle[2] = tempR[8];
      if (!std::isnan(euler_angle[0])) {
        iy = 0;
      } else {
        iy = -1;
        i = 2;
        exitg1 = false;
        while ((!exitg1) && (i < 4)) {
          if (!std::isnan(euler_angle[i - 1])) {
            iy = i - 1;
            exitg1 = true;
          } else {
            i++;
          }
        }
      }

      if (iy + 1 == 0) {
        x_size = 0;
      } else {
        scale = euler_angle[iy];
        x_size = iy;
        for (i = iy + 2; i < 4; i++) {
          if (scale < euler_angle[i - 1]) {
            scale = euler_angle[i - 1];
            x_size = i - 1;
          }
        }
      }

      scale = std::fmod(static_cast<real_T>(x_size) + 1.0, 3.0);
      absxk = std::fmod((static_cast<real_T>(x_size) + 1.0) + 1.0, 3.0);
      t = std::sqrt(((tempR[3 * x_size + x_size] - tempR[3 * static_cast<int32_T>
                      (scale) + static_cast<int32_T>(scale)]) - tempR[3 *
                     static_cast<int32_T>(absxk) + static_cast<int32_T>(absxk)])
                    + 1.0);
      euler_angle[0] = 0.0;
      euler_angle[1] = 0.0;
      euler_angle[2] = 0.0;
      euler_angle[x_size] = t / 2.0;
      euler_angle[static_cast<int32_T>(scale)] = (tempR[3 * x_size +
        static_cast<int32_T>(scale)] + tempR[3 * static_cast<int32_T>(scale) +
        x_size]) / (2.0 * t);
      euler_angle[static_cast<int32_T>(absxk)] = (tempR[3 * x_size +
        static_cast<int32_T>(absxk)] + tempR[3 * static_cast<int32_T>(absxk) +
        x_size]) / (2.0 * t);
      scale = 3.3121686421112381E-170;
      absxk = std::abs(euler_angle[0]);
      if (absxk > 3.3121686421112381E-170) {
        q_idx_3 = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        q_idx_3 = t * t;
      }

      absxk = std::abs(euler_angle[1]);
      if (absxk > scale) {
        t = scale / absxk;
        q_idx_3 = q_idx_3 * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        q_idx_3 += t * t;
      }

      absxk = std::abs(euler_angle[2]);
      if (absxk > scale) {
        t = scale / absxk;
        q_idx_3 = q_idx_3 * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        q_idx_3 += t * t;
      }

      q_idx_3 = scale * std::sqrt(q_idx_3);
      st_idx_1 = euler_angle[0];
      st_idx_1 = st_idx_2 * st_idx_1 / q_idx_3;
      euler_angle[0] = st_idx_1;
      st_idx_1 = euler_angle[1];
      st_idx_1 = st_idx_2 * st_idx_1 / q_idx_3;
      euler_angle[1] = st_idx_1;
    }

    uavNavigationRL2D_2022a_icra_B->roll_ref = euler_angle[0];
    uavNavigationRL2D_2022a_icra_B->pitch_ref = euler_angle[1];
    uavNavigationRL2D_2022a_icra_B->f = q_idx_1;

    /* End of MATLAB Function: '<S33>/attitude mapper ' */

    /* ZeroOrderHold: '<S4>/Zero-Order Hold5' incorporates:
     *  MATLAB Function: '<S10>/R_I_B constructor'
     */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[2] == 0) {
      /* ZeroOrderHold: '<S4>/Zero-Order Hold5' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold5 =
        uavNavigationRL2D_2022a_icra_B->TransportDelay6;

      /* ZeroOrderHold: '<S4>/Zero-Order Hold4' */
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold4 =
        uavNavigationRL2D_2022a_icra_B->TransportDelay11;

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      euler_angle[0] = uavNavigationRL2D_2022a_icra_B->ZeroOrderHold5;
      euler_angle[1] = uavNavigationRL2D_2022a_icra_B->ZeroOrderHold4;
      scale = 3.3121686421112381E-170;
      st_idx_1 = euler_angle[0];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      absxk = std::abs(st_idx_1);
      if (absxk > 3.3121686421112381E-170) {
        st_idx_2 = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        st_idx_2 = t * t;
      }

      euler_angle[0] = 0.0;
      st_idx_1 = euler_angle[1];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      absxk = std::abs(st_idx_1);
      if (absxk > scale) {
        t = scale / absxk;
        st_idx_2 = st_idx_2 * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        st_idx_2 += t * t;
      }

      euler_angle[1] = 0.0;
      euler_angle[2] = 0.0;

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      st_idx_2 = scale * std::sqrt(st_idx_2);
      if (st_idx_2 != 0.0) {
        euler_angle[0] = uavNavigationRL2D_2022a_icra_B->ZeroOrderHold5 /
          st_idx_2;
        euler_angle[1] = uavNavigationRL2D_2022a_icra_B->ZeroOrderHold4 /
          st_idx_2;
        euler_angle[2] = 0.0 / st_idx_2;
      }

      st_idx_2 = std::fmin(1.5707963267948966, std::fmax(st_idx_2,
        -1.5707963267948966));
      q_idx_1 = std::sin(st_idx_2 / 2.0);
      q_idx_2 = std::cos(st_idx_2 / 2.0);
      b_y[0] = q_idx_2;
      b_y[1] = euler_angle[0] * q_idx_1;
      b_y[2] = euler_angle[1] * q_idx_1;
      b_y[3] = euler_angle[2] * q_idx_1;
      q_idx_3 = b_y[0];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      st_idx_2 = q_idx_3 * q_idx_3;
      q_idx_3 = b_y[1];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      q_idx_1 = q_idx_3 * q_idx_3;
      q_idx_3 = b_y[2];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      q_idx_2 = q_idx_3 * q_idx_3;
      q_idx_3 = b_y[3];

      /* MATLAB Function: '<S10>/R_I_B constructor' */
      q_idx_3 *= q_idx_3;
      st_idx_2 += q_idx_1;
      st_idx_2 += q_idx_2;
      st_idx_2 += q_idx_3;
      scale = 1.0 / std::sqrt(st_idx_2);
      st_idx_2 = b_y[0] * scale;
      q_idx_1 = b_y[1] * scale;
      q_idx_2 = b_y[2] * scale;
      q_idx_3 = b_y[3] * scale;
      roll_pitch_DCM[0] = 1.0 - (q_idx_2 * q_idx_2 + q_idx_3 * q_idx_3) * 2.0;
      roll_pitch_DCM[1] = (q_idx_1 * q_idx_2 - st_idx_2 * q_idx_3) * 2.0;
      roll_pitch_DCM[2] = (q_idx_1 * q_idx_3 + st_idx_2 * q_idx_2) * 2.0;
      roll_pitch_DCM[3] = (q_idx_1 * q_idx_2 + st_idx_2 * q_idx_3) * 2.0;
      roll_pitch_DCM[4] = 1.0 - (q_idx_1 * q_idx_1 + q_idx_3 * q_idx_3) * 2.0;
      roll_pitch_DCM[5] = (q_idx_2 * q_idx_3 - st_idx_2 * q_idx_1) * 2.0;
      roll_pitch_DCM[6] = (q_idx_1 * q_idx_3 - st_idx_2 * q_idx_2) * 2.0;
      roll_pitch_DCM[7] = (q_idx_2 * q_idx_3 + st_idx_2 * q_idx_1) * 2.0;
      roll_pitch_DCM[8] = 1.0 - (q_idx_1 * q_idx_1 + q_idx_2 * q_idx_2) * 2.0;
      std::memcpy(&tempR[0], &roll_pitch_DCM[0], 9U * sizeof(real_T));
      q_idx_2 = std::cos(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_3 = std::sin(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      q_idx_1 = std::sin(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      st_idx_2 = std::cos(uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3);
      a[0] = q_idx_2;
      a[3] = -q_idx_3;
      a[6] = 0.0;
      a[1] = q_idx_1;
      a[4] = st_idx_2;
      a[7] = 0.0;
      for (iy = 0; iy < 3; iy++) {
        g_data = static_cast<int8_T>(iy + 1);
        roll_pitch_DCM[g_data - 1] = tempR[(g_data - 1) * 3];
        g_data = static_cast<int8_T>(iy + 1);
        roll_pitch_DCM[g_data + 2] = tempR[(g_data - 1) * 3 + 1];
        g_data = static_cast<int8_T>(iy + 1);
        roll_pitch_DCM[g_data + 5] = tempR[(g_data - 1) * 3 + 2];
        a[3 * iy + 2] = c[iy];
      }

      for (i = 0; i < 3; i++) {
        for (iy = 0; iy <= 0; iy += 2) {
          _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i],
                        _mm_set1_pd(0.0));
          tmp_0 = _mm_loadu_pd(&a[iy]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(roll_pitch_DCM[3 * i]), tmp_0);
          tmp_1 = _mm_loadu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + iy]);
          tmp_0 = _mm_add_pd(tmp_1, tmp_0);
          _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i],
                        tmp_0);
          tmp_0 = _mm_loadu_pd(&a[iy + 3]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(roll_pitch_DCM[3 * i + 1]), tmp_0);
          tmp_1 = _mm_loadu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + iy]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i],
                        tmp_0);
          tmp_0 = _mm_loadu_pd(&a[iy + 6]);
          tmp_0 = _mm_mul_pd(_mm_set1_pd(roll_pitch_DCM[3 * i + 2]), tmp_0);
          tmp_1 = _mm_loadu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + iy]);
          tmp_0 = _mm_add_pd(tmp_0, tmp_1);
          _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i],
                        tmp_0);
        }

        for (iy = 2; iy < 3; iy++) {
          uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i] = 0.0;
          uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i] =
            uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + iy] + roll_pitch_DCM[3
            * i] * a[iy];
          uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i] = roll_pitch_DCM[3 *
            i + 1] * a[iy + 3] + uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i +
            iy];
          uavNavigationRL2D_2022a_icra_B->DCM_p[iy + 3 * i] = roll_pitch_DCM[3 *
            i + 2] * a[iy + 6] + uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i +
            iy];
        }
      }
    }

    /* End of ZeroOrderHold: '<S4>/Zero-Order Hold5' */

    /* MATLAB Function: '<S33>/Body frame attitude error1' */
    for (i = 0; i < 3; i++) {
      for (iy = 0; iy < 3; iy++) {
        tempR[i + 3 * iy] = 0.0;
        tempR[i + 3 * iy] += uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i] *
          uavNavigationRL2D_2022a_icra_B->DCM_ref[3 * iy];
        tempR[i + 3 * iy] += uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + 1] *
          uavNavigationRL2D_2022a_icra_B->DCM_ref[3 * iy + 1];
        tempR[i + 3 * iy] += uavNavigationRL2D_2022a_icra_B->DCM_p[3 * i + 2] *
          uavNavigationRL2D_2022a_icra_B->DCM_ref[3 * iy + 2];
      }
    }

    t = tempR[0];
    t += tempR[4];
    t += tempR[8];
    st_idx_2 = std::acos(std::fmin(1.0, std::fmax((t - 1.0) / 2.0, -1.0)));
    if (st_idx_2 != 0.0) {
      uavNavigationRL2D_2022a_icra_B->e_roll = (tempR[5] - tempR[7]) * st_idx_2 /
        (2.0 * std::sin(st_idx_2));
      uavNavigationRL2D_2022a_icra_B->e_pitch = (tempR[6] - tempR[2]) * st_idx_2
        / (2.0 * std::sin(st_idx_2));
      uavNavigationRL2D_2022a_icra_B->e_yaw = (tempR[1] - tempR[3]) * st_idx_2 /
        (2.0 * std::sin(st_idx_2));
    } else {
      uavNavigationRL2D_2022a_icra_B->e_roll = 0.0;
      uavNavigationRL2D_2022a_icra_B->e_pitch = 0.0;
      uavNavigationRL2D_2022a_icra_B->e_yaw = 0.0;
    }

    /* End of MATLAB Function: '<S33>/Body frame attitude error1' */

    /* MATLAB Function: '<S4>/throttle limiter' incorporates:
     *  Constant: '<Root>/True3'
     */
    euler_angle[0] = uavNavigationRL2D_2022a_icra_B->e_roll;
    euler_angle[1] = uavNavigationRL2D_2022a_icra_B->e_pitch;
    scale = 3.3121686421112381E-170;
    absxk = std::abs(euler_angle[0]);
    if (absxk > 3.3121686421112381E-170) {
      st_idx_2 = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      st_idx_2 = t * t;
    }

    absxk = std::abs(euler_angle[1]);
    if (absxk > scale) {
      t = scale / absxk;
      st_idx_2 = st_idx_2 * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      st_idx_2 += t * t;
    }

    uavNavigationRL2D_2022a_icra_B->throttle_out = std::cos(scale * std::sqrt
      (st_idx_2)) * uavNavigationRL2D_2022a_icra_B->f;
    if ((uavNavigationRL2D_2022a_icra_P.True3_Value != 0.0) &&
        (uavNavigationRL2D_2022a_icra_B->throttle_out < 0.0)) {
      uavNavigationRL2D_2022a_icra_B->throttle_out = 0.0;
    }

    /* End of MATLAB Function: '<S4>/throttle limiter' */

    /* Product: '<S13>/Product2' incorporates:
     *  Constant: '<S13>/Constant3'
     */
    uavNavigationRL2D_2022a_icra_B->Product2_f = kp_roll *
      uavNavigationRL2D_2022a_icra_B->e_roll;

    /* TransportDelay: '<S5>/Transport Delay15' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay15_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_roll;
      uavNavigationRL2D_2022a_icra_B->TransportDelay15 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay15_InitOutput,
        0,
        0);
    }

    /* Gain: '<S13>/Gain3' */
    uavNavigationRL2D_2022a_icra_B->Gain3 =
      uavNavigationRL2D_2022a_icra_P.Gain3_Gain *
      uavNavigationRL2D_2022a_icra_B->TransportDelay15;

    /* Product: '<S13>/Product3' incorporates:
     *  Constant: '<S13>/Constant4'
     */
    uavNavigationRL2D_2022a_icra_B->Product3_d = kd_roll *
      uavNavigationRL2D_2022a_icra_B->Gain3;

    /* Sum: '<S13>/Sum1' */
    uavNavigationRL2D_2022a_icra_B->Sum1_f0 =
      uavNavigationRL2D_2022a_icra_B->Product2_f +
      uavNavigationRL2D_2022a_icra_B->Product3_d;

    /* Product: '<S12>/Product2' incorporates:
     *  Constant: '<S12>/Constant1'
     */
    uavNavigationRL2D_2022a_icra_B->Product2_nn = kp_pitch *
      uavNavigationRL2D_2022a_icra_B->e_pitch;

    /* TransportDelay: '<S5>/Transport Delay16' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay16_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_pitch;
      uavNavigationRL2D_2022a_icra_B->TransportDelay16 = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *uBuffer,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.CircularBufSize,
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Last,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail,
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head,
        uavNavigationRL2D_2022a_icra_P.TransportDelay16_InitOutput,
        0,
        0);
    }

    /* Gain: '<S12>/Gain3' */
    uavNavigationRL2D_2022a_icra_B->Gain3_g =
      uavNavigationRL2D_2022a_icra_P.Gain3_Gain_p *
      uavNavigationRL2D_2022a_icra_B->TransportDelay16;

    /* Product: '<S12>/Product3' incorporates:
     *  Constant: '<S12>/Constant2'
     */
    uavNavigationRL2D_2022a_icra_B->Product3_e = kd_pitch *
      uavNavigationRL2D_2022a_icra_B->Gain3_g;

    /* Sum: '<S12>/Sum1' */
    uavNavigationRL2D_2022a_icra_B->Sum1_a =
      uavNavigationRL2D_2022a_icra_B->Product2_nn +
      uavNavigationRL2D_2022a_icra_B->Product3_e;

    /* Gain: '<S14>/Gain2' */
    uavNavigationRL2D_2022a_icra_B->Gain2 = kd_yaw *
      uavNavigationRL2D_2022a_icra_B->Gain13;

    /* Gain: '<S14>/Gain3' */
    uavNavigationRL2D_2022a_icra_B->Gain3_b = kp_yaw *
      uavNavigationRL2D_2022a_icra_B->e_yaw;

    /* Sum: '<S14>/Sum6' */
    uavNavigationRL2D_2022a_icra_B->Sum6_h =
      uavNavigationRL2D_2022a_icra_B->Gain3_b -
      uavNavigationRL2D_2022a_icra_B->Gain2;

    /* Saturate: '<S4>/Maximum  Command Authority' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[0];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[0];
    if (uavNavigationRL2D_2022a_icra_B->throttle_out > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->throttle_out < q_idx_1)) {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->throttle_out;
    }

    /* Saturate: '<S4>/Maximum  Command Authority' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[0] = q_idx_1;

    /* Saturate: '<S4>/Maximum  Command Authority' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[1];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[1];
    if (uavNavigationRL2D_2022a_icra_B->Sum1_f0 > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->Sum1_f0 < q_idx_1)) {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->Sum1_f0;
    }

    /* Saturate: '<S4>/Maximum  Command Authority' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[1] = q_idx_1;

    /* Saturate: '<S4>/Maximum  Command Authority' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[2];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[2];
    if (uavNavigationRL2D_2022a_icra_B->Sum1_a > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->Sum1_a < q_idx_1)) {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->Sum1_a;
    }

    /* Saturate: '<S4>/Maximum  Command Authority' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[2] = q_idx_1;

    /* Saturate: '<S4>/Maximum  Command Authority' */
    q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[3];
    st_idx_2 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[3];
    if (uavNavigationRL2D_2022a_icra_B->Sum6_h > st_idx_2) {
      q_idx_1 = st_idx_2;
    } else if (!(uavNavigationRL2D_2022a_icra_B->Sum6_h < q_idx_1)) {
      q_idx_1 = uavNavigationRL2D_2022a_icra_B->Sum6_h;
    }

    /* Saturate: '<S4>/Maximum  Command Authority' */
    uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[3] = q_idx_1;

    /* Product: '<S4>/Product' incorporates:
     *  Constant: '<S4>/Constant5'
     */
    tmp_2 = &uavNavigationRL2D_2022a_icra_P.Constant5_Value[0];
    b_y[0] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[0];
    b_y[1] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[1];
    b_y[2] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[2];
    b_y[3] = uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[3];
    for (i = 0; i < 4; i++) {
      q_idx_1 = tmp_2[i] * b_y[0];
      q_idx_1 += tmp_2[i + 4] * b_y[1];
      q_idx_1 += tmp_2[i + 8] * b_y[2];
      q_idx_1 += tmp_2[i + 12] * b_y[3];

      /* Product: '<S4>/Product' */
      uavNavigationRL2D_2022a_icra_B->Product_p[i] = q_idx_1;

      /* Saturate: '<S4>/Maximum  Command Authority 1' */
      q_idx_3 = uavNavigationRL2D_2022a_icra_B->Product_p[i];
      q_idx_1 = uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_Lowe_f[i];
      st_idx_2 =
        uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_Uppe_c[i];
      if (q_idx_3 > st_idx_2) {
        q_idx_3 = st_idx_2;
      } else if (q_idx_3 < q_idx_1) {
        q_idx_3 = q_idx_1;
      }

      /* Saturate: '<S4>/Maximum  Command Authority 1' */
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[i] = q_idx_3;
    }

    /* End of Product: '<S4>/Product' */

    /* TransportDelay: '<S3>/Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK_j.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_act;
      if (tau_act == 0.0)
        uavNavigationRL2D_2022a_icra_B->TransportDelay_m =
          uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[0];
      else
        uavNavigationRL2D_2022a_icra_B->TransportDelay_m = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *uBuffer,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.CircularBufSize,
          &uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Last,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail,
          uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head,
          uavNavigationRL2D_2022a_icra_P.TransportDelay_InitOutput_g,
          0,
          0);
    }

    /* Sum: '<S3>/Sum' */
    uavNavigationRL2D_2022a_icra_B->Sum =
      uavNavigationRL2D_2022a_icra_B->TransportDelay_m -
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited;

    /* Gain: '<S3>/Gain' */
    st_idx_2 = 1.0 / T_motor;

    /* Gain: '<S3>/Gain' */
    uavNavigationRL2D_2022a_icra_B->Gain_m = st_idx_2 *
      uavNavigationRL2D_2022a_icra_B->Sum;

    /* TransportDelay: '<S3>/Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK_e.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_act;
      if (tau_act == 0.0)
        uavNavigationRL2D_2022a_icra_B->TransportDelay1_i =
          uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[1];
      else
        uavNavigationRL2D_2022a_icra_B->TransportDelay1_i = rt_TDelayInterpolate
          (
           tMinusDelay,
           0.0,
           *uBuffer,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.CircularBufSize,
           &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Last,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head,
           uavNavigationRL2D_2022a_icra_P.TransportDelay1_InitOutput_n,
           0,
           0);
    }

    /* Sum: '<S3>/Sum1' */
    uavNavigationRL2D_2022a_icra_B->Sum1_d =
      uavNavigationRL2D_2022a_icra_B->TransportDelay1_i -
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited1;

    /* Gain: '<S3>/Gain1' */
    st_idx_2 = 1.0 / T_motor;

    /* Gain: '<S3>/Gain1' */
    uavNavigationRL2D_2022a_icra_B->Gain1 = st_idx_2 *
      uavNavigationRL2D_2022a_icra_B->Sum1_d;

    /* TransportDelay: '<S3>/Transport Delay2' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK_c.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_act;
      if (tau_act == 0.0)
        uavNavigationRL2D_2022a_icra_B->TransportDelay2_f =
          uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[2];
      else
        uavNavigationRL2D_2022a_icra_B->TransportDelay2_f = rt_TDelayInterpolate
          (
           tMinusDelay,
           0.0,
           *uBuffer,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.CircularBufSize,
           &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Last,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head,
           uavNavigationRL2D_2022a_icra_P.TransportDelay2_InitOutput_i,
           0,
           0);
    }

    /* Sum: '<S3>/Sum2' */
    uavNavigationRL2D_2022a_icra_B->Sum2 =
      uavNavigationRL2D_2022a_icra_B->TransportDelay2_f -
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited2;

    /* Gain: '<S3>/Gain2' */
    st_idx_2 = 1.0 / T_motor;

    /* Gain: '<S3>/Gain2' */
    uavNavigationRL2D_2022a_icra_B->Gain2_a = st_idx_2 *
      uavNavigationRL2D_2022a_icra_B->Sum2;

    /* TransportDelay: '<S3>/Transport Delay3' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK_k.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      real_T tMinusDelay = simTime - tau_act;
      if (tau_act == 0.0)
        uavNavigationRL2D_2022a_icra_B->TransportDelay3_m =
          uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[3];
      else
        uavNavigationRL2D_2022a_icra_B->TransportDelay3_m = rt_TDelayInterpolate
          (
           tMinusDelay,
           0.0,
           *uBuffer,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.CircularBufSize,
           &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Last,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail,
           uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head,
           uavNavigationRL2D_2022a_icra_P.TransportDelay3_InitOutput_i,
           0,
           0);
    }

    /* Sum: '<S3>/Sum3' */
    uavNavigationRL2D_2022a_icra_B->Sum3_o =
      uavNavigationRL2D_2022a_icra_B->TransportDelay3_m -
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited3;

    /* Gain: '<S3>/Gain3' */
    st_idx_2 = 1.0 / T_motor;

    /* Gain: '<S3>/Gain3' */
    uavNavigationRL2D_2022a_icra_B->Gain3_d = st_idx_2 *
      uavNavigationRL2D_2022a_icra_B->Sum3_o;
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      for (i = 0; i < 3; i++) {
        /* Math: '<S40>/Transpose' */
        uavNavigationRL2D_2022a_icra_B->Transpose[3 * i] =
          uavNavigationRL2D_2022a_icra_B->DCM[i];
        uavNavigationRL2D_2022a_icra_B->Transpose[3 * i + 1] =
          uavNavigationRL2D_2022a_icra_B->DCM[i + 3];
        uavNavigationRL2D_2022a_icra_B->Transpose[3 * i + 2] =
          uavNavigationRL2D_2022a_icra_B->DCM[i + 6];
      }
    }

    /* Integrator: '<S40>/Integrator5' */
    uavNavigationRL2D_2022a_icra_B->Integrator5_l[0] =
      uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[0];
    uavNavigationRL2D_2022a_icra_B->Integrator5_l[1] =
      uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[1];
    uavNavigationRL2D_2022a_icra_B->Integrator5_l[2] =
      uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[2];

    /* Product: '<S40>/MatrixMultiply' incorporates:
     *  Math: '<S40>/Transpose'
     */
    std::memcpy(&tempR[0], &uavNavigationRL2D_2022a_icra_B->Transpose[0], 9U *
                sizeof(real_T));
    euler_angle[0] = uavNavigationRL2D_2022a_icra_B->Integrator5_l[0];
    euler_angle[1] = uavNavigationRL2D_2022a_icra_B->Integrator5_l[1];
    euler_angle[2] = uavNavigationRL2D_2022a_icra_B->Integrator5_l[2];
    for (i = 0; i <= 0; i += 2) {
      /* Product: '<S40>/MatrixMultiply' */
      tmp_0 = _mm_loadu_pd(&tempR[i]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(euler_angle[0]));
      tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&tempR[i + 3]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[1]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S40>/MatrixMultiply' */
      tmp_1 = _mm_loadu_pd(&tempR[i + 6]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[2]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S40>/MatrixMultiply' */
      _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->MatrixMultiply[i], tmp_0);

      /* Gain: '<S40>/Gain7' incorporates:
       *  Product: '<S40>/MatrixMultiply'
       */
      tmp_0 = _mm_loadu_pd(&uavNavigationRL2D_2022a_icra_P.Gain7_Gain[i]);

      /* Gain: '<S40>/Gain7' incorporates:
       *  Product: '<S40>/MatrixMultiply'
       */
      tmp_1 = _mm_loadu_pd(&uavNavigationRL2D_2022a_icra_B->MatrixMultiply[i]);
      tmp_0 = _mm_mul_pd(tmp_0, tmp_1);

      /* Gain: '<S40>/Gain7' incorporates:
       *  Product: '<S40>/MatrixMultiply'
       */
      _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->Gain7[i], tmp_0);
    }

    for (i = 2; i < 3; i++) {
      /* Product: '<S40>/MatrixMultiply' */
      q_idx_1 = tempR[i] * euler_angle[0];
      q_idx_1 += tempR[i + 3] * euler_angle[1];
      q_idx_1 += tempR[i + 6] * euler_angle[2];

      /* Product: '<S40>/MatrixMultiply' */
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply[i] = q_idx_1;

      /* Gain: '<S40>/Gain7' */
      uavNavigationRL2D_2022a_icra_B->Gain7[i] =
        uavNavigationRL2D_2022a_icra_P.Gain7_Gain[i] *
        uavNavigationRL2D_2022a_icra_B->MatrixMultiply[i];
    }

    /* Product: '<S40>/MatrixMultiply1' */
    std::memcpy(&tempR[0], &uavNavigationRL2D_2022a_icra_B->DCM[0], 9U * sizeof
                (real_T));
    euler_angle[0] = uavNavigationRL2D_2022a_icra_B->Gain7[0];
    euler_angle[1] = uavNavigationRL2D_2022a_icra_B->Gain7[1];
    euler_angle[2] = uavNavigationRL2D_2022a_icra_B->Gain7[2];
    for (i = 0; i <= 0; i += 2) {
      /* Product: '<S40>/MatrixMultiply1' */
      tmp_0 = _mm_loadu_pd(&tempR[i]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(euler_angle[0]));
      tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&tempR[i + 3]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[1]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S40>/MatrixMultiply1' */
      tmp_1 = _mm_loadu_pd(&tempR[i + 6]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[2]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S40>/MatrixMultiply1' */
      _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[i], tmp_0);
    }

    for (i = 2; i < 3; i++) {
      /* Product: '<S40>/MatrixMultiply1' */
      q_idx_1 = tempR[i] * euler_angle[0];
      q_idx_1 += tempR[i + 3] * euler_angle[1];
      q_idx_1 += tempR[i + 6] * euler_angle[2];

      /* Product: '<S40>/MatrixMultiply1' */
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[i] = q_idx_1;
    }

    /* Gain: '<S5>/Gain7' */
    st_idx_2 = K_z / T_aero;

    /* Gain: '<S5>/Gain7' */
    uavNavigationRL2D_2022a_icra_B->Gain7_b = st_idx_2 *
      uavNavigationRL2D_2022a_icra_B->Product1[0];

    /* SignalConversion generated from: '<S41>/MatrixMultiply' incorporates:
     *  Constant: '<S41>/Constant1'
     *  Constant: '<S41>/Constant2'
     */
    uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[0] =
      uavNavigationRL2D_2022a_icra_P.Constant2_Value_d;
    uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[1] =
      uavNavigationRL2D_2022a_icra_P.Constant1_Value_k;
    uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[2] =
      uavNavigationRL2D_2022a_icra_B->Gain7_b;

    /* Product: '<S41>/MatrixMultiply' */
    std::memcpy(&tempR[0], &uavNavigationRL2D_2022a_icra_B->DCM[0], 9U * sizeof
                (real_T));
    euler_angle[0] =
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[0];
    euler_angle[1] =
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[1];
    euler_angle[2] =
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[2];
    for (i = 0; i <= 0; i += 2) {
      /* Product: '<S41>/MatrixMultiply' */
      tmp_0 = _mm_loadu_pd(&tempR[i]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(euler_angle[0]));
      tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&tempR[i + 3]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[1]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S41>/MatrixMultiply' */
      tmp_1 = _mm_loadu_pd(&tempR[i + 6]);
      tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(euler_angle[2]));
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);

      /* Product: '<S41>/MatrixMultiply' */
      _mm_storeu_pd(&uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[i], tmp_0);
    }

    for (i = 2; i < 3; i++) {
      /* Product: '<S41>/MatrixMultiply' */
      q_idx_1 = tempR[i] * euler_angle[0];
      q_idx_1 += tempR[i + 3] * euler_angle[1];
      q_idx_1 += tempR[i + 6] * euler_angle[2];

      /* Product: '<S41>/MatrixMultiply' */
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[i] = q_idx_1;
    }

    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Product: '<S5>/Product' incorporates:
       *  Constant: '<Root>/True3'
       *  Constant: '<S5>/Constant3'
       */
      uavNavigationRL2D_2022a_icra_B->Product_o =
        uavNavigationRL2D_2022a_icra_P.True3_Value *
        uavNavigationRL2D_2022a_icra_P.Constant3_Value_b;
    }

    /* Sum: '<S5>/Sum' */
    uavNavigationRL2D_2022a_icra_B->Sum_b =
      uavNavigationRL2D_2022a_icra_B->Product_o +
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[2];

    /* Sum: '<S40>/Sum' */
    uavNavigationRL2D_2022a_icra_B->acc[0] =
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[0] -
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[0];
    uavNavigationRL2D_2022a_icra_B->acc[1] =
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[1] -
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[1];
    uavNavigationRL2D_2022a_icra_B->acc[2] =
      uavNavigationRL2D_2022a_icra_B->Sum_b -
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[2];

    /* Gain: '<S5>/Gain10' */
    uavNavigationRL2D_2022a_icra_B->Gain10 = K_pitch *
      uavNavigationRL2D_2022a_icra_B->Product1[2];

    /* Gain: '<S5>/Gain8' */
    uavNavigationRL2D_2022a_icra_B->Gain8 = K_roll *
      uavNavigationRL2D_2022a_icra_B->Product1[1];
  }

  if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M)) {
    real_T *lastU;

    /* Update for TransportDelay: '<S5>/Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail+
                     1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->Integrator5_l[0];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->Integrator5_l[1];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay2' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->Integrator5_l[2];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay3' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->acc[0];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay4' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay4_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->acc[1];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay5' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay5_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->acc[2];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay6' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay6_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head+1)
         : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.CircularBufSize
             -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail
                     +1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->roll_h;
    }

    /* Update for TransportDelay: '<S5>/Transport Delay11' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay11_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head+
                   1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->pitch;
    }

    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for Memory: '<S41>/Memory' */
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[0] =
        uavNavigationRL2D_2022a_icra_B->quat[0];
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[1] =
        uavNavigationRL2D_2022a_icra_B->quat[1];
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[2] =
        uavNavigationRL2D_2022a_icra_B->quat[2];
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[3] =
        uavNavigationRL2D_2022a_icra_B->quat[3];
    }

    /* Update for TransportDelay: '<S5>/Transport Delay14' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay14_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head+
                   1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->yaw;
    }

    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for Memory: '<S22>/Memory2' */
      uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput =
        uavNavigationRL2D_2022a_icra_B->RelationalOperator;

      /* Update for Memory: '<S23>/Memory' */
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_l =
        uavNavigationRL2D_2022a_icra_B->MinMax;

      /* Update for Memory: '<S24>/Memory1' */
      uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput =
        uavNavigationRL2D_2022a_icra_B->MinMax1;
    }

    /* Update for Derivative: '<S5>/Derivative3' */
    if (uavNavigationRL2D_2022a_icra_DW->TimeStampA == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampB == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA <
               uavNavigationRL2D_2022a_icra_DW->TimeStampB) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA;
    } else {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB;
    }

    *lastU = uavNavigationRL2D_2022a_icra_B->Sum1[0];

    /* End of Update for Derivative: '<S5>/Derivative3' */

    /* Update for Derivative: '<S5>/Derivative4' */
    if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_b == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA_b =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_i;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampB_l == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB_l =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_h;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_b <
               uavNavigationRL2D_2022a_icra_DW->TimeStampB_l) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA_b =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_i;
    } else {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB_l =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_h;
    }

    *lastU = uavNavigationRL2D_2022a_icra_B->Sum1[1];

    /* End of Update for Derivative: '<S5>/Derivative4' */
    if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M) &&
        uavNavigationRL2D_2022a_icra_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for Memory: '<S25>/Memory2' */
      uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput_j =
        uavNavigationRL2D_2022a_icra_B->RelationalOperator_b;

      /* Update for Memory: '<S26>/Memory' */
      uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_k =
        uavNavigationRL2D_2022a_icra_B->MinMax_n;

      /* Update for Memory: '<S27>/Memory1' */
      uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput_o =
        uavNavigationRL2D_2022a_icra_B->MinMax1_h;
    }

    /* Update for Derivative: '<S5>/Derivative5' */
    if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_h == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA_h =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_p;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampB_a == (rtInf)) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB_a =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_o;
    } else if (uavNavigationRL2D_2022a_icra_DW->TimeStampA_h <
               uavNavigationRL2D_2022a_icra_DW->TimeStampB_a) {
      uavNavigationRL2D_2022a_icra_DW->TimeStampA_h =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_p;
    } else {
      uavNavigationRL2D_2022a_icra_DW->TimeStampB_a =
        uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      lastU = &uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_o;
    }

    *lastU = uavNavigationRL2D_2022a_icra_B->Sum1[2];

    /* End of Update for Derivative: '<S5>/Derivative5' */

    /* Update for TransportDelay: '<S5>/Transport Delay15' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay15_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head+
                   1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->TransferFcn7;
    }

    /* Update for TransportDelay: '<S5>/Transport Delay16' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay16_PWORK.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head+
                   1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head] =
        uavNavigationRL2D_2022a_icra_B->TransferFcn8;
    }

    /* Update for TransportDelay: '<S3>/Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK_j.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head+
                   1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head] = simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head] =
        uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[0];
    }

    /* Update for TransportDelay: '<S3>/Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK_e.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head
                   +1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head] =
        simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head] =
        uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[1];
    }

    /* Update for TransportDelay: '<S3>/Transport Delay2' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK_c.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head
                   +1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head] =
        simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head] =
        uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[2];
    }

    /* Update for TransportDelay: '<S3>/Transport Delay3' */
    {
      real_T **uBuffer = (real_T**)
        &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK_k.TUbufferPtrs[0];
      real_T simTime = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
      uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head =
        ((uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head <
          (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.CircularBufSize
           -1)) ? (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head
                   +1) : 0);
      if (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head ==
          uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail =
          ((uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail <
            (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.CircularBufSize
             -1)) ?
           (uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail+1) : 0);
      }

      (*uBuffer +
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.CircularBufSize)
        [uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head] =
        simTime;
      (*uBuffer)[uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head] =
        uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[3];
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(uavNavigationRL2D_2022a_icra_M)) {
    rt_ertODEUpdateContinuousStates(uavNavigationRL2D_2022a_icra_M->solverInfo,
      uavNavigationRL2D_2022a_icra_M);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++uavNavigationRL2D_2022a_icra_M->Timing.clockTick0)) {
      ++uavNavigationRL2D_2022a_icra_M->Timing.clockTickH0;
    }

    uavNavigationRL2D_2022a_icra_M->Timing.t[0] = rtsiGetSolverStopTime
      (uavNavigationRL2D_2022a_icra_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.0005s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.0005, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      uavNavigationRL2D_2022a_icra_M->Timing.clockTick1++;
      if (!uavNavigationRL2D_2022a_icra_M->Timing.clockTick1) {
        uavNavigationRL2D_2022a_icra_M->Timing.clockTickH1++;
      }
    }

    rate_scheduler(uavNavigationRL2D_2022a_icra_M);
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void uavNavigationRL2D_2022a_icra_derivatives(RT_MODEL_uavNavigationRL2D_20_T *
  const uavNavigationRL2D_2022a_icra_M)
{
  B_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_B{
    uavNavigationRL2D_2022a_icra_M->blockIO };

  X_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_X{
    uavNavigationRL2D_2022a_icra_M->contStates };

  XDot_uavNavigationRL2D_2022a__T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_uavNavigationRL2D_2022a__T *)
             uavNavigationRL2D_2022a_icra_M->derivs);

  /* Derivatives for Integrator: '<S5>/Integrator7' */
  _rtXdot->Integrator7_CSTATE = uavNavigationRL2D_2022a_icra_B->TransportDelay;

  /* Derivatives for Integrator: '<S5>/Integrator6' */
  _rtXdot->Integrator6_CSTATE = uavNavigationRL2D_2022a_icra_B->TransportDelay1;

  /* Derivatives for Integrator: '<S5>/Integrator5' */
  _rtXdot->Integrator5_CSTATE = uavNavigationRL2D_2022a_icra_B->TransportDelay2;

  /* Derivatives for TransferFcn: '<S5>/Transfer Fcn7' */
  _rtXdot->TransferFcn7_CSTATE = uavNavigationRL2D_2022a_icra_P.TransferFcn7_A *
    uavNavigationRL2D_2022a_icra_X->TransferFcn7_CSTATE;
  _rtXdot->TransferFcn7_CSTATE += uavNavigationRL2D_2022a_icra_B->Gain8;

  /* Derivatives for TransferFcn: '<S5>/Transfer Fcn8' */
  _rtXdot->TransferFcn8_CSTATE = uavNavigationRL2D_2022a_icra_P.TransferFcn8_A *
    uavNavigationRL2D_2022a_icra_X->TransferFcn8_CSTATE;
  _rtXdot->TransferFcn8_CSTATE += uavNavigationRL2D_2022a_icra_B->Gain10;

  /* Derivatives for Integrator: '<S3>/Integrator Limited' */
  lsat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE <= min_u);
  usat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE >= max_u);
  if (((!lsat) && (!usat)) || (lsat && (uavNavigationRL2D_2022a_icra_B->Gain_m >
        0.0)) || (usat && (uavNavigationRL2D_2022a_icra_B->Gain_m < 0.0))) {
    _rtXdot->IntegratorLimited_CSTATE = uavNavigationRL2D_2022a_icra_B->Gain_m;
  } else {
    /* in saturation */
    _rtXdot->IntegratorLimited_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator Limited' */

  /* Derivatives for Integrator: '<S3>/Integrator Limited1' */
  lsat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE <= min_u);
  usat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE >= max_u);
  if (((!lsat) && (!usat)) || (lsat && (uavNavigationRL2D_2022a_icra_B->Gain1 >
        0.0)) || (usat && (uavNavigationRL2D_2022a_icra_B->Gain1 < 0.0))) {
    _rtXdot->IntegratorLimited1_CSTATE = uavNavigationRL2D_2022a_icra_B->Gain1;
  } else {
    /* in saturation */
    _rtXdot->IntegratorLimited1_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator Limited1' */

  /* Derivatives for Integrator: '<S3>/Integrator Limited2' */
  lsat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE <= min_u);
  usat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE >= max_u);
  if (((!lsat) && (!usat)) || (lsat && (uavNavigationRL2D_2022a_icra_B->Gain2_a >
        0.0)) || (usat && (uavNavigationRL2D_2022a_icra_B->Gain2_a < 0.0))) {
    _rtXdot->IntegratorLimited2_CSTATE = uavNavigationRL2D_2022a_icra_B->Gain2_a;
  } else {
    /* in saturation */
    _rtXdot->IntegratorLimited2_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator Limited2' */

  /* Derivatives for Integrator: '<S3>/Integrator Limited3' */
  lsat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE <= min_u);
  usat = (uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE >= max_u);
  if (((!lsat) && (!usat)) || (lsat && (uavNavigationRL2D_2022a_icra_B->Gain3_d >
        0.0)) || (usat && (uavNavigationRL2D_2022a_icra_B->Gain3_d < 0.0))) {
    _rtXdot->IntegratorLimited3_CSTATE = uavNavigationRL2D_2022a_icra_B->Gain3_d;
  } else {
    /* in saturation */
    _rtXdot->IntegratorLimited3_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator Limited3' */

  /* Derivatives for TransferFcn: '<S22>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE = uavNavigationRL2D_2022a_icra_P.TransferFcn1_A *
    uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += uavNavigationRL2D_2022a_icra_B->Gain5;

  /* Derivatives for TransferFcn: '<S25>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE_p =
    uavNavigationRL2D_2022a_icra_P.TransferFcn1_A_d *
    uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE_p;
  _rtXdot->TransferFcn1_CSTATE_p += uavNavigationRL2D_2022a_icra_B->Gain5_h;

  /* Derivatives for Integrator: '<S40>/Integrator5' */
  _rtXdot->Integrator5_CSTATE_h[0] = uavNavigationRL2D_2022a_icra_B->acc[0];
  _rtXdot->Integrator5_CSTATE_h[1] = uavNavigationRL2D_2022a_icra_B->acc[1];
  _rtXdot->Integrator5_CSTATE_h[2] = uavNavigationRL2D_2022a_icra_B->acc[2];
}

/* Model initialize function */
void uavNavigationRL2D_2022a_icra_initialize(RT_MODEL_uavNavigationRL2D_20_T *
  const uavNavigationRL2D_2022a_icra_M)
{
  B_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_B{
    uavNavigationRL2D_2022a_icra_M->blockIO };

  DW_uavNavigationRL2D_2022a_ic_T *uavNavigationRL2D_2022a_icra_DW{
    uavNavigationRL2D_2022a_icra_M->dwork };

  X_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_X{
    uavNavigationRL2D_2022a_icra_M->contStates };

  /* Start for TransportDelay: '<S5>/Transport Delay' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK.CircularBufSize = 1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay1' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay1_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay2' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay2_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay3' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay3_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay4' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay4_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay4_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay5' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay5_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay5_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay6' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay6_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay6_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay11' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay11_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay11_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay14' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay14_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay14_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for Constant: '<S2>/Constant1' */
  uavNavigationRL2D_2022a_icra_B->Constant1 =
    uavNavigationRL2D_2022a_icra_P.Constant1_Value_d;

  /* Start for TransportDelay: '<S5>/Transport Delay15' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay15_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay15_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S5>/Transport Delay16' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay16_RWORK.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_IWORK.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay16_InitOutput;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_PWORK.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S3>/Transport Delay' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK_d.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_IWORK_k.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay_InitOutput_g;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_PWORK_j.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S3>/Transport Delay1' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK_o.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_IWORK_e.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay1_InitOutput_n;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_PWORK_e.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S3>/Transport Delay2' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK_p.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_IWORK_h.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay2_InitOutput_i;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_PWORK_c.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* Start for TransportDelay: '<S3>/Transport Delay3' */
  {
    real_T *pBuffer =
      &uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK_i.TUbufferArea[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Tail = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Head = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.Last = 0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_IWORK_n.CircularBufSize =
      1024;
    pBuffer[0] = uavNavigationRL2D_2022a_icra_P.TransportDelay3_InitOutput_i;
    pBuffer[1024] = uavNavigationRL2D_2022a_icra_M->Timing.t[0];
    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_PWORK_k.TUbufferPtrs[0] =
      (void *) &pBuffer[0];
  }

  /* InitializeConditions for Integrator: '<S5>/Integrator7' */
  uavNavigationRL2D_2022a_icra_X->Integrator7_CSTATE =
    uavNavigationRL2D_2022a_icra_P.Integrator7_IC;

  /* InitializeConditions for Integrator: '<S5>/Integrator6' */
  uavNavigationRL2D_2022a_icra_X->Integrator6_CSTATE =
    uavNavigationRL2D_2022a_icra_P.Integrator6_IC;

  /* InitializeConditions for Integrator: '<S5>/Integrator5' */
  uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE =
    uavNavigationRL2D_2022a_icra_P.Integrator5_IC;

  /* InitializeConditions for TransferFcn: '<S5>/Transfer Fcn7' */
  uavNavigationRL2D_2022a_icra_X->TransferFcn7_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Transfer Fcn8' */
  uavNavigationRL2D_2022a_icra_X->TransferFcn8_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator Limited' */
  uavNavigationRL2D_2022a_icra_X->IntegratorLimited_CSTATE =
    uavNavigationRL2D_2022a_icra_P.IntegratorLimited_IC;

  /* InitializeConditions for Integrator: '<S3>/Integrator Limited1' */
  uavNavigationRL2D_2022a_icra_X->IntegratorLimited1_CSTATE =
    uavNavigationRL2D_2022a_icra_P.IntegratorLimited1_IC;

  /* InitializeConditions for Integrator: '<S3>/Integrator Limited2' */
  uavNavigationRL2D_2022a_icra_X->IntegratorLimited2_CSTATE =
    uavNavigationRL2D_2022a_icra_P.IntegratorLimited2_IC;

  /* InitializeConditions for Integrator: '<S3>/Integrator Limited3' */
  uavNavigationRL2D_2022a_icra_X->IntegratorLimited3_CSTATE =
    uavNavigationRL2D_2022a_icra_P.IntegratorLimited3_IC;

  /* InitializeConditions for Memory: '<S41>/Memory' */
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[0] =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition[0];
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[1] =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition[1];
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[2] =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition[2];
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[3] =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition[3];

  /* InitializeConditions for TransferFcn: '<S22>/Transfer Fcn1' */
  uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE = 0.0;

  /* InitializeConditions for Memory: '<S22>/Memory2' */
  uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput =
    uavNavigationRL2D_2022a_icra_P.Memory2_InitialCondition;

  /* InitializeConditions for Memory: '<S23>/Memory' */
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_l =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition_c;

  /* InitializeConditions for Memory: '<S24>/Memory1' */
  uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput =
    uavNavigationRL2D_2022a_icra_P.Memory1_InitialCondition;

  /* InitializeConditions for Derivative: '<S5>/Derivative3' */
  uavNavigationRL2D_2022a_icra_DW->TimeStampA = (rtInf);
  uavNavigationRL2D_2022a_icra_DW->TimeStampB = (rtInf);

  /* InitializeConditions for Derivative: '<S5>/Derivative4' */
  uavNavigationRL2D_2022a_icra_DW->TimeStampA_b = (rtInf);
  uavNavigationRL2D_2022a_icra_DW->TimeStampB_l = (rtInf);

  /* InitializeConditions for TransferFcn: '<S25>/Transfer Fcn1' */
  uavNavigationRL2D_2022a_icra_X->TransferFcn1_CSTATE_p = 0.0;

  /* InitializeConditions for Memory: '<S25>/Memory2' */
  uavNavigationRL2D_2022a_icra_DW->Memory2_PreviousInput_j =
    uavNavigationRL2D_2022a_icra_P.Memory2_InitialCondition_g;

  /* InitializeConditions for Memory: '<S26>/Memory' */
  uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_k =
    uavNavigationRL2D_2022a_icra_P.Memory_InitialCondition_f;

  /* InitializeConditions for Memory: '<S27>/Memory1' */
  uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput_o =
    uavNavigationRL2D_2022a_icra_P.Memory1_InitialCondition_h;

  /* InitializeConditions for Derivative: '<S5>/Derivative5' */
  uavNavigationRL2D_2022a_icra_DW->TimeStampA_h = (rtInf);
  uavNavigationRL2D_2022a_icra_DW->TimeStampB_a = (rtInf);

  /* InitializeConditions for Integrator: '<S40>/Integrator5' */
  uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[0] =
    uavNavigationRL2D_2022a_icra_P.Integrator5_IC_k;
  uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[1] =
    uavNavigationRL2D_2022a_icra_P.Integrator5_IC_k;
  uavNavigationRL2D_2022a_icra_X->Integrator5_CSTATE_h[2] =
    uavNavigationRL2D_2022a_icra_P.Integrator5_IC_k;

  /* SystemInitialize for MATLAB Function: '<S41>/kinematic differential equation' */
  uavNavigationRL2D_2022a_icra_DW->last_t = 0.0;
}

/* Model terminate function */
void uavNavigationRL2D_2022a_icra_terminate(RT_MODEL_uavNavigationRL2D_20_T
  * uavNavigationRL2D_2022a_icra_M)
{
  rt_FREE(uavNavigationRL2D_2022a_icra_M->solverInfo);

  /* model code */
  rt_FREE(uavNavigationRL2D_2022a_icra_M->blockIO);
  rt_FREE(uavNavigationRL2D_2022a_icra_M->contStates);
  rt_FREE(uavNavigationRL2D_2022a_icra_M->dwork);
  delete uavNavigationRL2D_2022a_icra_M;
}

/* Model data allocation function */
RT_MODEL_uavNavigationRL2D_20_T *uavNavigationRL2D_2022a_icra(void)
{
  RT_MODEL_uavNavigationRL2D_20_T *uavNavigationRL2D_2022a_icra_M;
  uavNavigationRL2D_2022a_icra_M = new RT_MODEL_uavNavigationRL2D_20_T();
  if (uavNavigationRL2D_2022a_icra_M == (nullptr)) {
    return (nullptr);
  }

  {
    /* Setup solver object */
    RTWSolverInfo *rt_SolverInfo{ (RTWSolverInfo *) malloc(sizeof(RTWSolverInfo))
    };

    rt_VALIDATE_MEMORY(uavNavigationRL2D_2022a_icra_M,rt_SolverInfo);
    uavNavigationRL2D_2022a_icra_M->solverInfo = (rt_SolverInfo);
    rtsiSetSimTimeStepPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
                          &uavNavigationRL2D_2022a_icra_M->Timing.simTimeStep);
    rtsiSetTPtr(uavNavigationRL2D_2022a_icra_M->solverInfo, &rtmGetTPtr
                (uavNavigationRL2D_2022a_icra_M));
    rtsiSetStepSizePtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
                       &uavNavigationRL2D_2022a_icra_M->Timing.stepSize0);
    rtsiSetdXPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
                 &uavNavigationRL2D_2022a_icra_M->derivs);
    rtsiSetContStatesPtr(uavNavigationRL2D_2022a_icra_M->solverInfo, (real_T **)
                         &uavNavigationRL2D_2022a_icra_M->contStates);
    rtsiSetNumContStatesPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
      &uavNavigationRL2D_2022a_icra_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
      &uavNavigationRL2D_2022a_icra_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr
      (uavNavigationRL2D_2022a_icra_M->solverInfo,
       &uavNavigationRL2D_2022a_icra_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
      &uavNavigationRL2D_2022a_icra_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
                          (&rtmGetErrorStatus(uavNavigationRL2D_2022a_icra_M)));
    rtsiSetRTModelPtr(uavNavigationRL2D_2022a_icra_M->solverInfo,
                      uavNavigationRL2D_2022a_icra_M);
  }

  rtsiSetSolverName(uavNavigationRL2D_2022a_icra_M->solverInfo,"ode8");

  /* block I/O */
  {
    B_uavNavigationRL2D_2022a_icr_T *b{ (B_uavNavigationRL2D_2022a_icr_T *)
      malloc(sizeof(B_uavNavigationRL2D_2022a_icr_T)) };

    rt_VALIDATE_MEMORY(uavNavigationRL2D_2022a_icra_M,b);
    uavNavigationRL2D_2022a_icra_M->blockIO = (b);
  }

  /* states (continuous) */
  {
    X_uavNavigationRL2D_2022a_icr_T *x{ (X_uavNavigationRL2D_2022a_icr_T *)
      malloc(sizeof(X_uavNavigationRL2D_2022a_icr_T)) };

    rt_VALIDATE_MEMORY(uavNavigationRL2D_2022a_icra_M,x);
    uavNavigationRL2D_2022a_icra_M->contStates = (x);
  }

  /* states (dwork) */
  {
    DW_uavNavigationRL2D_2022a_ic_T *dwork{
      static_cast<DW_uavNavigationRL2D_2022a_ic_T *>(malloc(sizeof
      (DW_uavNavigationRL2D_2022a_ic_T))) };

    rt_VALIDATE_MEMORY(uavNavigationRL2D_2022a_icra_M,dwork);
    uavNavigationRL2D_2022a_icra_M->dwork = (dwork);
  }

  {
    B_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_B{
      uavNavigationRL2D_2022a_icra_M->blockIO };

    DW_uavNavigationRL2D_2022a_ic_T *uavNavigationRL2D_2022a_icra_DW{
      uavNavigationRL2D_2022a_icra_M->dwork };

    X_uavNavigationRL2D_2022a_icr_T *uavNavigationRL2D_2022a_icra_X{
      uavNavigationRL2D_2022a_icra_M->contStates };

    /* initialize non-finites */
    rt_InitInfAndNaN(sizeof(real_T));

    /* non-finite (run-time) assignments */
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[0] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[1] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[2] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_UpperS[3] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[0] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[1] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[2] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority1_LowerS[3] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[0] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[1] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[2] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_UpperSa[3] = rtInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[0] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[1] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[2] =
      rtMinusInf;
    uavNavigationRL2D_2022a_icra_P.MaximumCommandAuthority_LowerSa[3] =
      rtMinusInf;
    rtsiSetSimTimeStep(uavNavigationRL2D_2022a_icra_M->solverInfo,
                       MAJOR_TIME_STEP);
    uavNavigationRL2D_2022a_icra_M->intgData.deltaY=
      uavNavigationRL2D_2022a_icra_M->OdeDeltaY;
    uavNavigationRL2D_2022a_icra_M->intgData.f[0] =
      uavNavigationRL2D_2022a_icra_M->odeF[0];
    uavNavigationRL2D_2022a_icra_M->intgData.f[1] =
      uavNavigationRL2D_2022a_icra_M->odeF[1];
    uavNavigationRL2D_2022a_icra_M->intgData.f[2] =
      uavNavigationRL2D_2022a_icra_M->odeF[2];
    uavNavigationRL2D_2022a_icra_M->intgData.f[3] =
      uavNavigationRL2D_2022a_icra_M->odeF[3];
    uavNavigationRL2D_2022a_icra_M->intgData.f[4] =
      uavNavigationRL2D_2022a_icra_M->odeF[4];
    uavNavigationRL2D_2022a_icra_M->intgData.f[5] =
      uavNavigationRL2D_2022a_icra_M->odeF[5];
    uavNavigationRL2D_2022a_icra_M->intgData.f[6] =
      uavNavigationRL2D_2022a_icra_M->odeF[6];
    uavNavigationRL2D_2022a_icra_M->intgData.f[7] =
      uavNavigationRL2D_2022a_icra_M->odeF[7];
    uavNavigationRL2D_2022a_icra_M->intgData.f[8] =
      uavNavigationRL2D_2022a_icra_M->odeF[8];
    uavNavigationRL2D_2022a_icra_M->intgData.f[9] =
      uavNavigationRL2D_2022a_icra_M->odeF[9];
    uavNavigationRL2D_2022a_icra_M->intgData.f[10] =
      uavNavigationRL2D_2022a_icra_M->odeF[10];
    uavNavigationRL2D_2022a_icra_M->intgData.f[11] =
      uavNavigationRL2D_2022a_icra_M->odeF[11];
    uavNavigationRL2D_2022a_icra_M->intgData.f[12] =
      uavNavigationRL2D_2022a_icra_M->odeF[12];
    uavNavigationRL2D_2022a_icra_M->intgData.x0 =
      uavNavigationRL2D_2022a_icra_M->odeX0;
    uavNavigationRL2D_2022a_icra_M->contStates =
      ((X_uavNavigationRL2D_2022a_icr_T *) uavNavigationRL2D_2022a_icra_X);
    rtsiSetSolverData(uavNavigationRL2D_2022a_icra_M->solverInfo, static_cast<
                      void *>(&uavNavigationRL2D_2022a_icra_M->intgData));
    rtsiSetIsMinorTimeStepWithModeChange
      (uavNavigationRL2D_2022a_icra_M->solverInfo, false);
    rtmSetTPtr(uavNavigationRL2D_2022a_icra_M,
               &uavNavigationRL2D_2022a_icra_M->Timing.tArray[0]);
    uavNavigationRL2D_2022a_icra_M->Timing.stepSize0 = 0.0005;

    /* block I/O */
    (void) std::memset((static_cast<void *>(uavNavigationRL2D_2022a_icra_B)), 0,
                       sizeof(B_uavNavigationRL2D_2022a_icr_T));

    {
      int32_T i;
      for (i = 0; i < 16; i++) {
        uavNavigationRL2D_2022a_icra_B->MatrixDivide[i] = 0.0;
      }

      for (i = 0; i < 9; i++) {
        uavNavigationRL2D_2022a_icra_B->Transpose[i] = 0.0;
      }

      for (i = 0; i < 9; i++) {
        uavNavigationRL2D_2022a_icra_B->DCM[i] = 0.0;
      }

      for (i = 0; i < 9; i++) {
        uavNavigationRL2D_2022a_icra_B->DCM_ref[i] = 0.0;
      }

      for (i = 0; i < 9; i++) {
        uavNavigationRL2D_2022a_icra_B->DCM_p[i] = 0.0;
      }

      uavNavigationRL2D_2022a_icra_B->Integrator7 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Integrator6 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Integrator5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay4 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay6 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain = 0.0;
      uavNavigationRL2D_2022a_icra_B->roll = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay11 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain_n = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransferFcn7 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransferFcn8 = 0.0;
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited = 0.0;
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->IntegratorLimited3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product1[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product1[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product1[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product1[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain13 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Clock = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay14 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain_i = 0.0;
      uavNavigationRL2D_2022a_icra_B->Constant1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold8 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold9 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransferFcn1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory_k = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain4 = 0.0;
      uavNavigationRL2D_2022a_icra_B->MinMax = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->MinMax1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch1_e = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold2[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Derivative3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold10 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Derivative4 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold11 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum6 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_f = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch_b = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product1_o = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2_p = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum4 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3_o = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_k = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold6 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Error = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2_a = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain5_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransferFcn1_a = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory_g = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch_l = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain4_c = 0.0;
      uavNavigationRL2D_2022a_icra_B->MinMax_n = 0.0;
      uavNavigationRL2D_2022a_icra_B->Memory1_e = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch1_j = 0.0;
      uavNavigationRL2D_2022a_icra_B->MinMax1_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch1_g = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2_n = 0.0;
      uavNavigationRL2D_2022a_icra_B->Derivative5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold7 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Error_i = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3_p = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3_ox = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_hp = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold5 = 0.0;
      uavNavigationRL2D_2022a_icra_B->ZeroOrderHold4 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2_f = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay15 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain3 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3_d = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_f0 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product2_nn = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay16 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain3_g = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product3_e = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_a = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain3_b = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum6_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product_p[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product_p[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product_p[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product_p[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MaximumCommandAuthority1_m[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay_m = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain_m = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay1_i = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum1_d = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain1 = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay2_f = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain2_a = 0.0;
      uavNavigationRL2D_2022a_icra_B->TransportDelay3_m = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum3_o = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain3_d = 0.0;
      uavNavigationRL2D_2022a_icra_B->Integrator5_l[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Integrator5_l[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Integrator5_l[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain7[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain7[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain7[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply1[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain7_b = 0.0;
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->TmpSignalConversionAtMatrixMult[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->MatrixMultiply_i[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Product_o = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum_b = 0.0;
      uavNavigationRL2D_2022a_icra_B->acc[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->acc[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->acc[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain10 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Gain8 = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat_rate[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat_rate[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat_rate[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat_rate[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat[0] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat[1] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat[2] = 0.0;
      uavNavigationRL2D_2022a_icra_B->quat[3] = 0.0;
      uavNavigationRL2D_2022a_icra_B->roll_h = 0.0;
      uavNavigationRL2D_2022a_icra_B->pitch = 0.0;
      uavNavigationRL2D_2022a_icra_B->yaw = 0.0;
      uavNavigationRL2D_2022a_icra_B->throttle_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->f = 0.0;
      uavNavigationRL2D_2022a_icra_B->roll_ref = 0.0;
      uavNavigationRL2D_2022a_icra_B->pitch_ref = 0.0;
      uavNavigationRL2D_2022a_icra_B->e_roll = 0.0;
      uavNavigationRL2D_2022a_icra_B->e_pitch = 0.0;
      uavNavigationRL2D_2022a_icra_B->e_yaw = 0.0;
      uavNavigationRL2D_2022a_icra_B->x_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->y_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch2 = 0.0;
      uavNavigationRL2D_2022a_icra_B->Beta = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum4_l = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sign = 0.0;
      uavNavigationRL2D_2022a_icra_B->h = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch_p = 0.0;
      uavNavigationRL2D_2022a_icra_B->Switch2_f = 0.0;
      uavNavigationRL2D_2022a_icra_B->Beta_i = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sum4_p = 0.0;
      uavNavigationRL2D_2022a_icra_B->Sign_b = 0.0;
      uavNavigationRL2D_2022a_icra_B->h_d = 0.0;
      uavNavigationRL2D_2022a_icra_B->vel_cmd_x_trim = 0.0;
      uavNavigationRL2D_2022a_icra_B->vel_cmd_y_trim = 0.0;
      uavNavigationRL2D_2022a_icra_B->vel_cmd_z_trim = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.x_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_ld.y_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_l.x_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_l.y_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_n.x_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_n.y_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_c.x_out = 0.0;
      uavNavigationRL2D_2022a_icra_B->sf_referencecorrection_c.y_out = 0.0;
    }

    /* Storage classes */
    {
      int32_T i;
      for (i = 0; i < 12; i++) {
        uavNavigationRL2D_states_output[i] = 0.0;
      }
    }

    uavNavigati_velocity_references[0] = 0.0;
    uavNavigati_velocity_references[1] = 0.0;
    uavNavigati_velocity_references[2] = 0.0;

    /* Storage classes */
    uavNavigationRL2_velocity_cmd_x = 0.0;
    uavNavigationRL2_velocity_cmd_y = 0.0;
    uavNavigationRL2_velocity_cmd_z = 0.0;
    uavNavigationRL2D_2022a_u_cmd_z = 0.0;

    /* states (continuous) */
    {
      (void) std::memset(static_cast<void *>(uavNavigationRL2D_2022a_icra_X), 0,
                         sizeof(X_uavNavigationRL2D_2022a_icr_T));
    }

    /* states (dwork) */
    (void) std::memset(static_cast<void *>(uavNavigationRL2D_2022a_icra_DW), 0,
                       sizeof(DW_uavNavigationRL2D_2022a_ic_T));

    {
      int32_T i;
      for (i = 0; i < 16; i++) {
        uavNavigationRL2D_2022a_icra_DW->MatrixDivide_DWORK1[i] = 0.0;
      }
    }

    {
      int32_T i;
      for (i = 0; i < 16; i++) {
        uavNavigationRL2D_2022a_icra_DW->MatrixDivide_DWORK3[i] = 0.0;
      }
    }

    {
      int32_T i;
      for (i = 0; i < 16; i++) {
        uavNavigationRL2D_2022a_icra_DW->MatrixDivide_DWORK4[i] = 0.0;
      }
    }

    {
      int32_T i;
      for (i = 0; i < 16; i++) {
        uavNavigationRL2D_2022a_icra_DW->MatrixDivide_DWORK5[i] = 0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[0] = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[1] = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[2] = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput[3] = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_l = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampA = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampB = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampA_b = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_i = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampB_l = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_h = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory_PreviousInput_k = 0.0;
    uavNavigationRL2D_2022a_icra_DW->Memory1_PreviousInput_o = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampA_h = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeA_p = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TimeStampB_a = 0.0;
    uavNavigationRL2D_2022a_icra_DW->LastUAtTimeB_o = 0.0;
    uavNavigationRL2D_2022a_icra_DW->last_t = 0.0;
    uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay4_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay4_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay5_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay5_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay6_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay6_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay11_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay11_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay14_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay14_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay15_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay15_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay16_RWORK.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay16_RWORK.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK_d.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay_RWORK_d.TUbufferArea[i] =
          0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK_o.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay1_RWORK_o.TUbufferArea[i]
          = 0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK_p.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay2_RWORK_p.TUbufferArea[i]
          = 0.0;
      }
    }

    uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK_i.modelTStart = 0.0;

    {
      int32_T i;
      for (i = 0; i < 2048; i++) {
        uavNavigationRL2D_2022a_icra_DW->TransportDelay3_RWORK_i.TUbufferArea[i]
          = 0.0;
      }
    }
  }

  return uavNavigationRL2D_2022a_icra_M;
}
