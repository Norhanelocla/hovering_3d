/*
 * uavNavigationRL2D_2022a_icra.h
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

#ifndef RTW_HEADER_uavNavigationRL2D_2022a_icra_h_
#define RTW_HEADER_uavNavigationRL2D_2022a_icra_h_
#include <stdlib.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "uavNavigationRL2D_2022a_icra_types.h"

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rt_nonfinite.h"

}

#include <cfloat>
#include <cmath>
#include <cstring>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlockIO
#define rtmGetBlockIO(rtm)             ((rtm)->blockIO)
#endif

#ifndef rtmSetBlockIO
#define rtmSetBlockIO(rtm, val)        ((rtm)->blockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeDeltaY
#define rtmGetOdeDeltaY(rtm)           ((rtm)->OdeDeltaY)
#endif

#ifndef rtmSetOdeDeltaY
#define rtmSetOdeDeltaY(rtm, val)      ((rtm)->OdeDeltaY = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeX0
#define rtmGetOdeX0(rtm)               ((rtm)->odeX0)
#endif

#ifndef rtmSetOdeX0
#define rtmSetOdeX0(rtm, val)          ((rtm)->odeX0 = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRootDWork
#define rtmGetRootDWork(rtm)           ((rtm)->dwork)
#endif

#ifndef rtmSetRootDWork
#define rtmSetRootDWork(rtm, val)      ((rtm)->dwork = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#define uavNavigationRL2D_2022a_icra_M_TYPE RT_MODEL_uavNavigationRL2D_20_T

/* Block signals for system '<S19>/reference correction' */
struct B_referencecorrection_uavNavi_T {
  real_T x_out;                        /* '<S19>/reference correction' */
  real_T y_out;                        /* '<S19>/reference correction' */
};

/* Block signals (default storage) */
struct B_uavNavigationRL2D_2022a_icr_T {
  real_T Integrator7;                  /* '<S5>/Integrator7' */
  real_T Integrator6;                  /* '<S5>/Integrator6' */
  real_T Integrator5;                  /* '<S5>/Integrator5' */
  real_T Sum1[3];                      /* '<S5>/Sum1' */
  real_T TransportDelay;               /* '<S5>/Transport Delay' */
  real_T TransportDelay1;              /* '<S5>/Transport Delay1' */
  real_T TransportDelay2;              /* '<S5>/Transport Delay2' */
  real_T TransportDelay3;              /* '<S5>/Transport Delay3' */
  real_T TransportDelay4;              /* '<S5>/Transport Delay4' */
  real_T TransportDelay5;              /* '<S5>/Transport Delay5' */
  real_T TransportDelay6;              /* '<S5>/Transport Delay6' */
  real_T Gain;                         /* '<S6>/Gain' */
  real_T roll;                         /* '<S2>/Gain' */
  real_T TransportDelay11;             /* '<S5>/Transport Delay11' */
  real_T Gain_n;                       /* '<S7>/Gain' */
  real_T TransferFcn7;                 /* '<S5>/Transfer Fcn7' */
  real_T TransferFcn8;                 /* '<S5>/Transfer Fcn8' */
  real_T MatrixDivide[16];             /* '<S5>/Matrix Divide' */
  real_T IntegratorLimited;            /* '<S3>/Integrator Limited' */
  real_T IntegratorLimited1;           /* '<S3>/Integrator Limited1' */
  real_T IntegratorLimited2;           /* '<S3>/Integrator Limited2' */
  real_T IntegratorLimited3;           /* '<S3>/Integrator Limited3' */
  real_T MaximumCommandAuthority1[4];  /* '<S5>/Maximum  Command Authority1' */
  real_T Product1[4];                  /* '<S5>/Product1' */
  real_T Gain13;                       /* '<S5>/Gain13' */
  real_T Clock;                        /* '<S41>/Clock' */
  real_T Memory[4];                    /* '<S41>/Memory' */
  real_T TransportDelay14;             /* '<S5>/Transport Delay14' */
  real_T ZeroOrderHold3;               /* '<S4>/Zero-Order Hold3' */
  real_T Gain_i;                       /* '<S10>/Gain' */
  real_T Constant1;                    /* '<S2>/Constant1' */
  real_T ZeroOrderHold8;               /* '<S4>/Zero-Order Hold8' */
  real_T ZeroOrderHold9;               /* '<S4>/Zero-Order Hold9' */
  real_T Sum3;                         /* '<S9>/Sum3' */
  real_T Product;                      /* '<S9>/Product' */
  real_T Gain5;                        /* '<S22>/Gain5' */
  real_T TransferFcn1;                 /* '<S22>/Transfer Fcn1' */
  real_T Memory_k;                     /* '<S23>/Memory' */
  real_T Switch;                       /* '<S23>/Switch' */
  real_T Gain4;                        /* '<S22>/Gain4' */
  real_T MinMax;                       /* '<S23>/MinMax' */
  real_T Memory1;                      /* '<S24>/Memory1' */
  real_T Switch1;                      /* '<S24>/Switch1' */
  real_T MinMax1;                      /* '<S24>/MinMax1' */
  real_T Switch1_e;                    /* '<S15>/Switch1' */
  real_T Product2;                     /* '<S15>/Product2' */
  real_T ZeroOrderHold2[3];            /* '<S4>/Zero-Order Hold2' */
  real_T Derivative3;                  /* '<S5>/Derivative3' */
  real_T ZeroOrderHold10;              /* '<S4>/Zero-Order Hold10' */
  real_T Derivative4;                  /* '<S5>/Derivative4' */
  real_T ZeroOrderHold11;              /* '<S4>/Zero-Order Hold11' */
  real_T Sum6;                         /* '<S9>/Sum6' */
  real_T Product3;                     /* '<S15>/Product3' */
  real_T Sum1_f;                       /* '<S15>/Sum1' */
  real_T Switch_b;                     /* '<S15>/Switch' */
  real_T Sum5;                         /* '<S9>/Sum5' */
  real_T Product1_o;                   /* '<S9>/Product1' */
  real_T Product2_p;                   /* '<S16>/Product2' */
  real_T Sum4;                         /* '<S9>/Sum4' */
  real_T Product3_o;                   /* '<S16>/Product3' */
  real_T Sum1_k;                       /* '<S16>/Sum1' */
  real_T ZeroOrderHold6;               /* '<S4>/Zero-Order Hold6' */
  real_T Error;                        /* '<S9>/Sum9' */
  real_T Product2_a;                   /* '<S9>/Product2' */
  real_T Gain5_h;                      /* '<S25>/Gain5' */
  real_T TransferFcn1_a;               /* '<S25>/Transfer Fcn1' */
  real_T Memory_g;                     /* '<S26>/Memory' */
  real_T Switch_l;                     /* '<S26>/Switch' */
  real_T Gain4_c;                      /* '<S25>/Gain4' */
  real_T MinMax_n;                     /* '<S26>/MinMax' */
  real_T Memory1_e;                    /* '<S27>/Memory1' */
  real_T Switch1_j;                    /* '<S27>/Switch1' */
  real_T MinMax1_h;                    /* '<S27>/MinMax1' */
  real_T Switch1_g;                    /* '<S17>/Switch1' */
  real_T Product2_n;                   /* '<S17>/Product2' */
  real_T Derivative5;                  /* '<S5>/Derivative5' */
  real_T ZeroOrderHold7;               /* '<S4>/Zero-Order Hold7' */
  real_T Error_i;                      /* '<S9>/Sum2' */
  real_T Product3_p;                   /* '<S17>/Product3' */
  real_T Sum1_h;                       /* '<S17>/Sum1' */
  real_T Switch_h;                     /* '<S9>/Switch' */
  real_T Product3_ox;                  /* '<S9>/Product3' */
  real_T Sum1_hp;                      /* '<S9>/Sum1' */
  real_T ZeroOrderHold5;               /* '<S4>/Zero-Order Hold5' */
  real_T ZeroOrderHold4;               /* '<S4>/Zero-Order Hold4' */
  real_T Product2_f;                   /* '<S13>/Product2' */
  real_T TransportDelay15;             /* '<S5>/Transport Delay15' */
  real_T Gain3;                        /* '<S13>/Gain3' */
  real_T Product3_d;                   /* '<S13>/Product3' */
  real_T Sum1_f0;                      /* '<S13>/Sum1' */
  real_T Product2_nn;                  /* '<S12>/Product2' */
  real_T TransportDelay16;             /* '<S5>/Transport Delay16' */
  real_T Gain3_g;                      /* '<S12>/Gain3' */
  real_T Product3_e;                   /* '<S12>/Product3' */
  real_T Sum1_a;                       /* '<S12>/Sum1' */
  real_T Gain2;                        /* '<S14>/Gain2' */
  real_T Gain3_b;                      /* '<S14>/Gain3' */
  real_T Sum6_h;                       /* '<S14>/Sum6' */
  real_T MaximumCommandAuthority[4];   /* '<S4>/Maximum  Command Authority' */
  real_T Product_p[4];                 /* '<S4>/Product' */
  real_T MaximumCommandAuthority1_m[4];/* '<S4>/Maximum  Command Authority 1' */
  real_T TransportDelay_m;             /* '<S3>/Transport Delay' */
  real_T Sum;                          /* '<S3>/Sum' */
  real_T Gain_m;                       /* '<S3>/Gain' */
  real_T TransportDelay1_i;            /* '<S3>/Transport Delay1' */
  real_T Sum1_d;                       /* '<S3>/Sum1' */
  real_T Gain1;                        /* '<S3>/Gain1' */
  real_T TransportDelay2_f;            /* '<S3>/Transport Delay2' */
  real_T Sum2;                         /* '<S3>/Sum2' */
  real_T Gain2_a;                      /* '<S3>/Gain2' */
  real_T TransportDelay3_m;            /* '<S3>/Transport Delay3' */
  real_T Sum3_o;                       /* '<S3>/Sum3' */
  real_T Gain3_d;                      /* '<S3>/Gain3' */
  real_T Transpose[9];                 /* '<S40>/Transpose' */
  real_T Integrator5_l[3];             /* '<S40>/Integrator5' */
  real_T MatrixMultiply[3];            /* '<S40>/MatrixMultiply' */
  real_T Gain7[3];                     /* '<S40>/Gain7' */
  real_T MatrixMultiply1[3];           /* '<S40>/MatrixMultiply1' */
  real_T Gain7_b;                      /* '<S5>/Gain7' */
  real_T TmpSignalConversionAtMatrixMult[3];
  real_T MatrixMultiply_i[3];          /* '<S41>/MatrixMultiply' */
  real_T Product_o;                    /* '<S5>/Product' */
  real_T Sum_b;                        /* '<S5>/Sum' */
  real_T acc[3];                       /* '<S40>/Sum' */
  real_T Gain10;                       /* '<S5>/Gain10' */
  real_T Gain8;                        /* '<S5>/Gain8' */
  real_T quat_rate[4];             /* '<S41>/kinematic differential equation' */
  real_T quat[4];                  /* '<S41>/kinematic differential equation' */
  real_T DCM[9];                   /* '<S41>/kinematic differential equation' */
  real_T roll_h;                   /* '<S41>/kinematic differential equation' */
  real_T pitch;                    /* '<S41>/kinematic differential equation' */
  real_T yaw;                      /* '<S41>/kinematic differential equation' */
  real_T throttle_out;                 /* '<S4>/throttle limiter' */
  real_T f;                            /* '<S33>/attitude mapper ' */
  real_T DCM_ref[9];                   /* '<S33>/attitude mapper ' */
  real_T roll_ref;                     /* '<S33>/attitude mapper ' */
  real_T pitch_ref;                    /* '<S33>/attitude mapper ' */
  real_T e_roll;                       /* '<S33>/Body frame attitude error1' */
  real_T e_pitch;                      /* '<S33>/Body frame attitude error1' */
  real_T e_yaw;                        /* '<S33>/Body frame attitude error1' */
  real_T DCM_p[9];                     /* '<S10>/R_I_B constructor' */
  real_T x_out;                        /* '<S18>/reference correction' */
  real_T y_out;                        /* '<S18>/reference correction' */
  real_T Switch2;                      /* '<S22>/Switch2' */
  real_T Beta;                         /* '<S22>/Beta' */
  real_T Sum4_l;                       /* '<S22>/Sum4' */
  real_T Sign;                         /* '<S22>/Sign' */
  real_T h;                            /* '<S22>/h' */
  real_T Switch_p;                     /* '<S17>/Switch' */
  real_T Switch2_f;                    /* '<S25>/Switch2' */
  real_T Beta_i;                       /* '<S25>/Beta' */
  real_T Sum4_p;                       /* '<S25>/Sum4' */
  real_T Sign_b;                       /* '<S25>/Sign' */
  real_T h_d;                          /* '<S25>/h' */
  real_T vel_cmd_x_trim;               /* '<Root>/MATLAB Function' */
  real_T vel_cmd_y_trim;               /* '<Root>/MATLAB Function' */
  real_T vel_cmd_z_trim;               /* '<Root>/MATLAB Function' */
  boolean_T RelationalOperator;        /* '<S22>/Relational Operator' */
  boolean_T Memory2;                   /* '<S22>/Memory2' */
  boolean_T LogicalOperator1;          /* '<S22>/Logical Operator1' */
  boolean_T LogicalOperator;           /* '<S22>/Logical Operator' */
  boolean_T LogicalOperator4;          /* '<S22>/Logical Operator4' */
  boolean_T LogicalOperator3;          /* '<S22>/Logical Operator3' */
  boolean_T RelationalOperator_b;      /* '<S25>/Relational Operator' */
  boolean_T Memory2_o;                 /* '<S25>/Memory2' */
  boolean_T LogicalOperator1_m;        /* '<S25>/Logical Operator1' */
  boolean_T LogicalOperator_c;         /* '<S25>/Logical Operator' */
  boolean_T LogicalOperator4_l;        /* '<S25>/Logical Operator4' */
  boolean_T LogicalOperator3_f;        /* '<S25>/Logical Operator3' */
  B_referencecorrection_uavNavi_T sf_referencecorrection_ld;/* '<S34>/reference correction' */
  B_referencecorrection_uavNavi_T sf_referencecorrection_l;/* '<S21>/reference correction' */
  B_referencecorrection_uavNavi_T sf_referencecorrection_n;/* '<S20>/reference correction' */
  B_referencecorrection_uavNavi_T sf_referencecorrection_c;/* '<S19>/reference correction' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_uavNavigationRL2D_2022a_ic_T {
  real_T MatrixDivide_DWORK1[16];      /* '<S5>/Matrix Divide' */
  real_T MatrixDivide_DWORK3[16];      /* '<S5>/Matrix Divide' */
  real_T MatrixDivide_DWORK4[16];      /* '<S5>/Matrix Divide' */
  real_T MatrixDivide_DWORK5[16];      /* '<S5>/Matrix Divide' */
  real_T Memory_PreviousInput[4];      /* '<S41>/Memory' */
  real_T Memory_PreviousInput_l;       /* '<S23>/Memory' */
  real_T Memory1_PreviousInput;        /* '<S24>/Memory1' */
  real_T TimeStampA;                   /* '<S5>/Derivative3' */
  real_T LastUAtTimeA;                 /* '<S5>/Derivative3' */
  real_T TimeStampB;                   /* '<S5>/Derivative3' */
  real_T LastUAtTimeB;                 /* '<S5>/Derivative3' */
  real_T TimeStampA_b;                 /* '<S5>/Derivative4' */
  real_T LastUAtTimeA_i;               /* '<S5>/Derivative4' */
  real_T TimeStampB_l;                 /* '<S5>/Derivative4' */
  real_T LastUAtTimeB_h;               /* '<S5>/Derivative4' */
  real_T Memory_PreviousInput_k;       /* '<S26>/Memory' */
  real_T Memory1_PreviousInput_o;      /* '<S27>/Memory1' */
  real_T TimeStampA_h;                 /* '<S5>/Derivative5' */
  real_T LastUAtTimeA_p;               /* '<S5>/Derivative5' */
  real_T TimeStampB_a;                 /* '<S5>/Derivative5' */
  real_T LastUAtTimeB_o;               /* '<S5>/Derivative5' */
  real_T last_t;                   /* '<S41>/kinematic differential equation' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK;              /* '<S5>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay1_RWORK;             /* '<S5>/Transport Delay1' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay2_RWORK;             /* '<S5>/Transport Delay2' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay3_RWORK;             /* '<S5>/Transport Delay3' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay4_RWORK;             /* '<S5>/Transport Delay4' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay5_RWORK;             /* '<S5>/Transport Delay5' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay6_RWORK;             /* '<S5>/Transport Delay6' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay11_RWORK;            /* '<S5>/Transport Delay11' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay14_RWORK;            /* '<S5>/Transport Delay14' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay15_RWORK;            /* '<S5>/Transport Delay15' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay16_RWORK;            /* '<S5>/Transport Delay16' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_d;            /* '<S3>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay1_RWORK_o;           /* '<S3>/Transport Delay1' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay2_RWORK_p;           /* '<S3>/Transport Delay2' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay3_RWORK_i;           /* '<S3>/Transport Delay3' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S5>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay1_PWORK;             /* '<S5>/Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay2_PWORK;             /* '<S5>/Transport Delay2' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay3_PWORK;             /* '<S5>/Transport Delay3' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay4_PWORK;             /* '<S5>/Transport Delay4' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay5_PWORK;             /* '<S5>/Transport Delay5' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay6_PWORK;             /* '<S5>/Transport Delay6' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay11_PWORK;            /* '<S5>/Transport Delay11' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay14_PWORK;            /* '<S5>/Transport Delay14' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay15_PWORK;            /* '<S5>/Transport Delay15' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay16_PWORK;            /* '<S5>/Transport Delay16' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_j;            /* '<S3>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay1_PWORK_e;           /* '<S3>/Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay2_PWORK_c;           /* '<S3>/Transport Delay2' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay3_PWORK_k;           /* '<S3>/Transport Delay3' */

  int32_T MatrixDivide_DWORK2[4];      /* '<S5>/Matrix Divide' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S5>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay1_IWORK;             /* '<S5>/Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay2_IWORK;             /* '<S5>/Transport Delay2' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay3_IWORK;             /* '<S5>/Transport Delay3' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay4_IWORK;             /* '<S5>/Transport Delay4' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay5_IWORK;             /* '<S5>/Transport Delay5' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay6_IWORK;             /* '<S5>/Transport Delay6' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay11_IWORK;            /* '<S5>/Transport Delay11' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay14_IWORK;            /* '<S5>/Transport Delay14' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay15_IWORK;            /* '<S5>/Transport Delay15' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay16_IWORK;            /* '<S5>/Transport Delay16' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_k;            /* '<S3>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay1_IWORK_e;           /* '<S3>/Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay2_IWORK_h;           /* '<S3>/Transport Delay2' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay3_IWORK_n;           /* '<S3>/Transport Delay3' */

  boolean_T Memory2_PreviousInput;     /* '<S22>/Memory2' */
  boolean_T Memory2_PreviousInput_j;   /* '<S25>/Memory2' */
};

/* Continuous states (default storage) */
struct X_uavNavigationRL2D_2022a_icr_T {
  real_T Integrator7_CSTATE;           /* '<S5>/Integrator7' */
  real_T Integrator6_CSTATE;           /* '<S5>/Integrator6' */
  real_T Integrator5_CSTATE;           /* '<S5>/Integrator5' */
  real_T TransferFcn7_CSTATE;          /* '<S5>/Transfer Fcn7' */
  real_T TransferFcn8_CSTATE;          /* '<S5>/Transfer Fcn8' */
  real_T IntegratorLimited_CSTATE;     /* '<S3>/Integrator Limited' */
  real_T IntegratorLimited1_CSTATE;    /* '<S3>/Integrator Limited1' */
  real_T IntegratorLimited2_CSTATE;    /* '<S3>/Integrator Limited2' */
  real_T IntegratorLimited3_CSTATE;    /* '<S3>/Integrator Limited3' */
  real_T TransferFcn1_CSTATE;          /* '<S22>/Transfer Fcn1' */
  real_T TransferFcn1_CSTATE_p;        /* '<S25>/Transfer Fcn1' */
  real_T Integrator5_CSTATE_h[3];      /* '<S40>/Integrator5' */
};

/* State derivatives (default storage) */
struct XDot_uavNavigationRL2D_2022a__T {
  real_T Integrator7_CSTATE;           /* '<S5>/Integrator7' */
  real_T Integrator6_CSTATE;           /* '<S5>/Integrator6' */
  real_T Integrator5_CSTATE;           /* '<S5>/Integrator5' */
  real_T TransferFcn7_CSTATE;          /* '<S5>/Transfer Fcn7' */
  real_T TransferFcn8_CSTATE;          /* '<S5>/Transfer Fcn8' */
  real_T IntegratorLimited_CSTATE;     /* '<S3>/Integrator Limited' */
  real_T IntegratorLimited1_CSTATE;    /* '<S3>/Integrator Limited1' */
  real_T IntegratorLimited2_CSTATE;    /* '<S3>/Integrator Limited2' */
  real_T IntegratorLimited3_CSTATE;    /* '<S3>/Integrator Limited3' */
  real_T TransferFcn1_CSTATE;          /* '<S22>/Transfer Fcn1' */
  real_T TransferFcn1_CSTATE_p;        /* '<S25>/Transfer Fcn1' */
  real_T Integrator5_CSTATE_h[3];      /* '<S40>/Integrator5' */
};

/* State disabled  */
struct XDis_uavNavigationRL2D_2022a__T {
  boolean_T Integrator7_CSTATE;        /* '<S5>/Integrator7' */
  boolean_T Integrator6_CSTATE;        /* '<S5>/Integrator6' */
  boolean_T Integrator5_CSTATE;        /* '<S5>/Integrator5' */
  boolean_T TransferFcn7_CSTATE;       /* '<S5>/Transfer Fcn7' */
  boolean_T TransferFcn8_CSTATE;       /* '<S5>/Transfer Fcn8' */
  boolean_T IntegratorLimited_CSTATE;  /* '<S3>/Integrator Limited' */
  boolean_T IntegratorLimited1_CSTATE; /* '<S3>/Integrator Limited1' */
  boolean_T IntegratorLimited2_CSTATE; /* '<S3>/Integrator Limited2' */
  boolean_T IntegratorLimited3_CSTATE; /* '<S3>/Integrator Limited3' */
  boolean_T TransferFcn1_CSTATE;       /* '<S22>/Transfer Fcn1' */
  boolean_T TransferFcn1_CSTATE_p;     /* '<S25>/Transfer Fcn1' */
  boolean_T Integrator5_CSTATE_h[3];   /* '<S40>/Integrator5' */
};

#ifndef ODE8_INTG
#define ODE8_INTG

/* ODE8 Integration Data */
struct ODE8_IntgData {
  real_T *deltaY;                      /* output diff */
  real_T *f[13];                       /* derivatives */
  real_T *x0;                          /* Initial State */
};

#endif

/* Parameters (default storage) */
struct P_uavNavigationRL2D_2022a_icr_T_ {
  real_T Beta_Gain;                    /* Expression: -0.73
                                        * Referenced by: '<S25>/Beta'
                                        */
  real_T h_Gain;                       /* Expression: 3
                                        * Referenced by: '<S25>/h'
                                        */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<S17>/Switch'
                                        */
  real_T Constant_Value;               /* Expression: 1
                                        * Referenced by: '<S15>/Constant'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S15>/Constant1'
                                        */
  real_T Beta_Gain_m;                  /* Expression: -0.73
                                        * Referenced by: '<S22>/Beta'
                                        */
  real_T h_Gain_k;                     /* Expression: 8
                                        * Referenced by: '<S22>/h'
                                        */
  real_T Constant_Value_m;             /* Expression: 1
                                        * Referenced by: '<S17>/Constant'
                                        */
  real_T Constant1_Value_j;            /* Expression: 0
                                        * Referenced by: '<S17>/Constant1'
                                        */
  real_T Integrator7_IC;               /* Expression: 0
                                        * Referenced by: '<S5>/Integrator7'
                                        */
  real_T Integrator6_IC;               /* Expression: 0
                                        * Referenced by: '<S5>/Integrator6'
                                        */
  real_T Integrator5_IC;               /* Expression: 0
                                        * Referenced by: '<S5>/Integrator5'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay'
                                        */
  real_T TransportDelay1_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay1'
                                        */
  real_T TransportDelay2_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay2'
                                        */
  real_T TransportDelay3_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay3'
                                        */
  real_T TransportDelay4_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay4'
                                        */
  real_T TransportDelay5_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay5'
                                        */
  real_T TransportDelay6_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay6'
                                        */
  real_T Gain_Gain;                    /* Expression: 180/pi
                                        * Referenced by: '<S6>/Gain'
                                        */
  real_T Gain_Gain_h;                  /* Expression: -1
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T TransportDelay11_InitOutput;  /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay11'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 180/pi
                                        * Referenced by: '<S7>/Gain'
                                        */
  real_T TransferFcn7_A;               /* Computed Parameter: TransferFcn7_A
                                        * Referenced by: '<S5>/Transfer Fcn7'
                                        */
  real_T TransferFcn7_C;               /* Computed Parameter: TransferFcn7_C
                                        * Referenced by: '<S5>/Transfer Fcn7'
                                        */
  real_T TransferFcn8_A;               /* Computed Parameter: TransferFcn8_A
                                        * Referenced by: '<S5>/Transfer Fcn8'
                                        */
  real_T TransferFcn8_C;               /* Computed Parameter: TransferFcn8_C
                                        * Referenced by: '<S5>/Transfer Fcn8'
                                        */
  real_T Constant2_Value[16];          /* Expression: eye(4)
                                        * Referenced by: '<S5>/Constant2'
                                        */
  real_T Constant4_Value[16];
                          /* Expression: [ 0.0489   -0.4581    0.5566    5.0787;
                             0.0489   -0.4581   -0.5566   -5.0787;
                             0.0489    0.4581    0.5566   -5.0787;
                             0.0489    0.4581   -0.5566    5.0787]
                           * Referenced by: '<S5>/Constant4'
                           */
  real_T IntegratorLimited_IC;         /* Expression: 0
                                        * Referenced by: '<S3>/Integrator Limited'
                                        */
  real_T IntegratorLimited1_IC;        /* Expression: 0
                                        * Referenced by: '<S3>/Integrator Limited1'
                                        */
  real_T IntegratorLimited2_IC;        /* Expression: 0
                                        * Referenced by: '<S3>/Integrator Limited2'
                                        */
  real_T IntegratorLimited3_IC;        /* Expression: 0
                                        * Referenced by: '<S3>/Integrator Limited3'
                                        */
  real_T MaximumCommandAuthority1_UpperS[4];/* Expression: [Inf Inf Inf Inf]
                                             * Referenced by: '<S5>/Maximum  Command Authority1'
                                             */
  real_T MaximumCommandAuthority1_LowerS[4];/* Expression: [-Inf -Inf -Inf -Inf]
                                             * Referenced by: '<S5>/Maximum  Command Authority1'
                                             */
  real_T Gain13_Gain;                  /* Expression: 1
                                        * Referenced by: '<S5>/Gain13'
                                        */
  real_T Memory_InitialCondition[4];   /* Expression: [1 0 0 0]
                                        * Referenced by: '<S41>/Memory'
                                        */
  real_T TransportDelay14_Delay;       /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay14'
                                        */
  real_T TransportDelay14_InitOutput;  /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay14'
                                        */
  real_T upper_limit_x_Value;          /* Expression: 1
                                        * Referenced by: '<Root>/upper_limit_x'
                                        */
  real_T upper_limit_y_Value;          /* Expression: 1
                                        * Referenced by: '<Root>/upper_limit_y'
                                        */
  real_T upper_limit_z_Value;          /* Expression: 5
                                        * Referenced by: '<Root>/upper_limit_z'
                                        */
  real_T lower_limit_x_Value;          /* Expression: -1
                                        * Referenced by: '<Root>/lower_limit_x'
                                        */
  real_T lower_limit_y_Value;          /* Expression: -1
                                        * Referenced by: '<Root>/lower_limit_y'
                                        */
  real_T lower_limit_z_Value;          /* Expression: 0
                                        * Referenced by: '<Root>/lower_limit_z'
                                        */
  real_T Constant5_Value[16];
                          /* Expression: [ 0.0489   -0.4581    0.5566    5.0787;
                             0.0489   -0.4581   -0.5566   -5.0787;
                             0.0489    0.4581    0.5566   -5.0787;
                             0.0489    0.4581   -0.5566    5.0787]
                           * Referenced by: '<S4>/Constant5'
                           */
  real_T Gain_Gain_f;                  /* Expression: -1
                                        * Referenced by: '<S10>/Gain'
                                        */
  real_T True1_Value;                  /* Expression: 1
                                        * Referenced by: '<Root>/True1'
                                        */
  real_T Constant1_Value_d;            /* Expression: 0
                                        * Referenced by: '<S2>/Constant1'
                                        */
  real_T Constant2_Value_e;            /* Expression: 0
                                        * Referenced by: '<S23>/Constant2'
                                        */
  real_T Gain5_Gain;                   /* Expression: -1
                                        * Referenced by: '<S22>/Gain5'
                                        */
  real_T TransferFcn1_A;               /* Computed Parameter: TransferFcn1_A
                                        * Referenced by: '<S22>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C;               /* Computed Parameter: TransferFcn1_C
                                        * Referenced by: '<S22>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_D;               /* Computed Parameter: TransferFcn1_D
                                        * Referenced by: '<S22>/Transfer Fcn1'
                                        */
  real_T Constant4_Value_a;            /* Expression: 0
                                        * Referenced by: '<S22>/Constant4'
                                        */
  real_T Memory_InitialCondition_c;    /* Expression: 0
                                        * Referenced by: '<S23>/Memory'
                                        */
  real_T Gain4_Gain;                   /* Expression: 1
                                        * Referenced by: '<S22>/Gain4'
                                        */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<S24>/Constant3'
                                        */
  real_T Memory1_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S24>/Memory1'
                                        */
  real_T Switch1_Threshold;            /* Expression: 5
                                        * Referenced by: '<S15>/Switch1'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 0
                                        * Referenced by: '<S15>/Switch'
                                        */
  real_T True2_Value;                  /* Expression: 1
                                        * Referenced by: '<Root>/True2'
                                        */
  real_T False1_Value;                 /* Expression: 0
                                        * Referenced by: '<Root>/False1'
                                        */
  real_T Constant2_Value_f;            /* Expression: 0
                                        * Referenced by: '<S26>/Constant2'
                                        */
  real_T Gain5_Gain_e;                 /* Expression: -1
                                        * Referenced by: '<S25>/Gain5'
                                        */
  real_T TransferFcn1_A_d;             /* Computed Parameter: TransferFcn1_A_d
                                        * Referenced by: '<S25>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C_b;             /* Computed Parameter: TransferFcn1_C_b
                                        * Referenced by: '<S25>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_D_d;             /* Computed Parameter: TransferFcn1_D_d
                                        * Referenced by: '<S25>/Transfer Fcn1'
                                        */
  real_T Constant4_Value_f;            /* Expression: 0
                                        * Referenced by: '<S25>/Constant4'
                                        */
  real_T Memory_InitialCondition_f;    /* Expression: 0
                                        * Referenced by: '<S26>/Memory'
                                        */
  real_T Gain4_Gain_c;                 /* Expression: 1
                                        * Referenced by: '<S25>/Gain4'
                                        */
  real_T Constant3_Value_a;            /* Expression: 0
                                        * Referenced by: '<S27>/Constant3'
                                        */
  real_T Memory1_InitialCondition_h;   /* Expression: 0
                                        * Referenced by: '<S27>/Memory1'
                                        */
  real_T Switch1_Threshold_n;          /* Expression: 5
                                        * Referenced by: '<S17>/Switch1'
                                        */
  real_T False2_Value;                 /* Expression: 0
                                        * Referenced by: '<Root>/False2'
                                        */
  real_T Switch_Threshold_m;           /* Expression: 0.5
                                        * Referenced by: '<S9>/Switch'
                                        */
  real_T True3_Value;                  /* Expression: 1
                                        * Referenced by: '<Root>/True3'
                                        */
  real_T max_tilt_ref_limit_deg_Value; /* Expression: 15
                                        * Referenced by: '<S2>/max_tilt_ref_limit_deg'
                                        */
  real_T virtualantigravityjusttomakeatt;/* Expression: 9.81
                                          * Referenced by: '<S33>/virtual anti-gravity (just to make attitude defined)'
                                          */
  real_T TransportDelay15_InitOutput;  /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay15'
                                        */
  real_T Gain3_Gain;                   /* Expression: -1
                                        * Referenced by: '<S13>/Gain3'
                                        */
  real_T TransportDelay16_InitOutput;  /* Expression: 0
                                        * Referenced by: '<S5>/Transport Delay16'
                                        */
  real_T Gain3_Gain_p;                 /* Expression: -1
                                        * Referenced by: '<S12>/Gain3'
                                        */
  real_T MaximumCommandAuthority_UpperSa[4];/* Expression: [Inf Inf Inf Inf]
                                             * Referenced by: '<S4>/Maximum  Command Authority'
                                             */
  real_T MaximumCommandAuthority_LowerSa[4];/* Expression: [-inf -Inf -Inf -Inf]
                                             * Referenced by: '<S4>/Maximum  Command Authority'
                                             */
  real_T MaximumCommandAuthority1_Uppe_c[4];/* Expression: [max_u max_u max_u max_u]
                                             * Referenced by: '<S4>/Maximum  Command Authority 1'
                                             */
  real_T MaximumCommandAuthority1_Lowe_f[4];/* Expression: [min_u min_u min_u min_u]
                                             * Referenced by: '<S4>/Maximum  Command Authority 1'
                                             */
  real_T TransportDelay_InitOutput_g;  /* Expression: 0
                                        * Referenced by: '<S3>/Transport Delay'
                                        */
  real_T TransportDelay1_InitOutput_n; /* Expression: 0
                                        * Referenced by: '<S3>/Transport Delay1'
                                        */
  real_T TransportDelay2_InitOutput_i; /* Expression: 0
                                        * Referenced by: '<S3>/Transport Delay2'
                                        */
  real_T TransportDelay3_InitOutput_i; /* Expression: 0
                                        * Referenced by: '<S3>/Transport Delay3'
                                        */
  real_T Integrator5_IC_k;             /* Expression: 0
                                        * Referenced by: '<S40>/Integrator5'
                                        */
  real_T Gain7_Gain[3];                /* Expression: [0, 0, 1/T_aero]
                                        * Referenced by: '<S40>/Gain7'
                                        */
  real_T Constant2_Value_d;            /* Expression: 0
                                        * Referenced by: '<S41>/Constant2'
                                        */
  real_T Constant1_Value_k;            /* Expression: 0
                                        * Referenced by: '<S41>/Constant1'
                                        */
  real_T Constant3_Value_b;            /* Expression: -9.81
                                        * Referenced by: '<S5>/Constant3'
                                        */
  boolean_T Memory2_InitialCondition;
                                 /* Computed Parameter: Memory2_InitialCondition
                                  * Referenced by: '<S22>/Memory2'
                                  */
  boolean_T Memory2_InitialCondition_g;
                               /* Computed Parameter: Memory2_InitialCondition_g
                                * Referenced by: '<S25>/Memory2'
                                */
};

/* Real-time Model Data Structure */
struct tag_RTM_uavNavigationRL2D_202_T {
  const char_T *errorStatus;
  RTWSolverInfo *solverInfo;
  B_uavNavigationRL2D_2022a_icr_T *blockIO;
  X_uavNavigationRL2D_2022a_icr_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_uavNavigationRL2D_2022a__T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T OdeDeltaY[14];
  real_T odeF[13][14];
  real_T odeX0[14];
  ODE8_IntgData intgData;
  DW_uavNavigationRL2D_2022a_ic_T *dwork;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block parameters (default storage) */
#ifdef __cplusplus

extern "C"
{

#endif

  extern P_uavNavigationRL2D_2022a_icr_T uavNavigationRL2D_2022a_icra_P;

#ifdef __cplusplus

}

#endif

/* External data declarations for dependent source files */
#ifdef __cplusplus

extern "C"
{

#endif

  extern const char *RT_MEMORY_ALLOCATION_ERROR;

#ifdef __cplusplus

}

#endif

extern P_uavNavigationRL2D_2022a_icr_T uavNavigationRL2D_2022a_icra_P;/* parameters */

#ifdef __cplusplus

extern "C"
{

#endif

  /* Model entry point functions */
  extern RT_MODEL_uavNavigationRL2D_20_T *uavNavigationRL2D_2022a_icra(void);
  extern void uavNavigationRL2D_2022a_icra_initialize
    (RT_MODEL_uavNavigationRL2D_20_T *const uavNavigationRL2D_2022a_icra_M);
  extern void uavNavigationRL2D_2022a_icra_step(RT_MODEL_uavNavigationRL2D_20_T *
    const uavNavigationRL2D_2022a_icra_M);
  extern void uavNavigationRL2D_2022a_icra_terminate
    (RT_MODEL_uavNavigationRL2D_20_T * uavNavigationRL2D_2022a_icra_M);

#ifdef __cplusplus

}

#endif

/* Exported data declaration */

/* Data with Exported storage */
extern real_T K_pitch;                 /* Referenced by: '<S5>/Gain10' */
extern real_T K_roll;                  /* Referenced by: '<S5>/Gain8' */
extern real_T K_z;                     /* Referenced by:
                                        * '<S5>/Constant'
                                        * '<S5>/Gain7'
                                        */
extern real_T T_aero;                  /* Referenced by:
                                        * '<S5>/Constant'
                                        * '<S5>/Gain7'
                                        */
extern real_T T_motor;                 /* Referenced by:
                                        * '<S3>/Gain'
                                        * '<S3>/Gain1'
                                        * '<S3>/Gain2'
                                        * '<S3>/Gain3'
                                        */
extern real_T kd_pitch;                /* Referenced by: '<S12>/Constant2' */
extern real_T kd_roll;                 /* Referenced by: '<S13>/Constant4' */
extern real_T kd_x;                    /* Referenced by: '<S15>/Constant5' */
extern real_T kd_y;                    /* Referenced by: '<S16>/Constant4' */
extern real_T kd_yaw;                  /* Referenced by: '<S14>/Gain2' */
extern real_T kd_z;                    /* Referenced by: '<S17>/Constant4' */
extern real_T kp_pitch;                /* Referenced by: '<S12>/Constant1' */
extern real_T kp_roll;                 /* Referenced by: '<S13>/Constant3' */
extern real_T kp_x;                    /* Referenced by: '<S15>/Constant2' */
extern real_T kp_y;                    /* Referenced by: '<S16>/Constant3' */
extern real_T kp_yaw;                  /* Referenced by: '<S14>/Gain3' */
extern real_T kp_z;                    /* Referenced by: '<S17>/Constant3' */
extern real_T max_u;                   /* Referenced by:
                                        * '<S3>/Integrator Limited'
                                        * '<S3>/Integrator Limited1'
                                        * '<S3>/Integrator Limited2'
                                        * '<S3>/Integrator Limited3'
                                        */
extern real_T min_u;                   /* Referenced by:
                                        * '<S3>/Integrator Limited'
                                        * '<S3>/Integrator Limited1'
                                        * '<S3>/Integrator Limited2'
                                        * '<S3>/Integrator Limited3'
                                        */
extern real_T tau_act;                 /* Referenced by:
                                        * '<S3>/Transport Delay'
                                        * '<S3>/Transport Delay1'
                                        * '<S3>/Transport Delay2'
                                        * '<S3>/Transport Delay3'
                                        */
extern real_T tau_imu;                 /* Referenced by:
                                        * '<S5>/Transport Delay3'
                                        * '<S5>/Transport Delay4'
                                        * '<S5>/Transport Delay5'
                                        */
extern real_T tau_pitch;               /* Referenced by:
                                        * '<S5>/Transport Delay11'
                                        * '<S5>/Transport Delay16'
                                        */
extern real_T tau_roll;                /* Referenced by:
                                        * '<S5>/Transport Delay15'
                                        * '<S5>/Transport Delay6'
                                        */
extern real_T tau_x;                 /* Referenced by: '<S5>/Transport Delay' */
extern real_T tau_y;                /* Referenced by: '<S5>/Transport Delay1' */
extern real_T tau_z;                /* Referenced by: '<S5>/Transport Delay2' */
extern real_T uavNavigati_velocity_references[3];/* '<Root>/velocity_references' */
extern real_T uavNavigationRL2D_2022a_u_cmd_z;/* '<Root>/u_cmd_z' */
extern real_T uavNavigationRL2D_states_output[12];/* '<Root>/states_output' */
extern real_T uavNavigationRL2_velocity_cmd_x;/* '<Root>/velocity_cmd_x' */
extern real_T uavNavigationRL2_velocity_cmd_y;/* '<Root>/velocity_cmd_y' */
extern real_T uavNavigationRL2_velocity_cmd_z;/* '<Root>/velocity_cmd_z' */
extern real_T uav_initial_x_position;  /* Referenced by: '<S5>/Constant1' */
extern real_T uav_initial_y_position;  /* Referenced by: '<S5>/Constant5' */
extern real_T uav_initial_z_position;  /* Referenced by: '<S5>/Constant6' */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'uavNavigationRL2D_2022a_icra'
 * '<S1>'   : 'uavNavigationRL2D_2022a_icra/MATLAB Function'
 * '<S2>'   : 'uavNavigationRL2D_2022a_icra/uavController2'
 * '<S3>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Actuator Dynamics1'
 * '<S4>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Controller'
 * '<S5>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Quadcopter Dynamics1'
 * '<S6>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Radians to Degrees'
 * '<S7>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Radians to Degrees1'
 * '<S8>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Inner loop Control'
 * '<S9>'   : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control '
 * '<S10>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization'
 * '<S11>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/throttle limiter'
 * '<S12>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Inner loop Control/pitch_controller'
 * '<S13>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Inner loop Control/roll controller '
 * '<S14>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Inner loop Control/yaw controller'
 * '<S15>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller x'
 * '<S16>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller y'
 * '<S17>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller z'
 * '<S18>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 1'
 * '<S19>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 2'
 * '<S20>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 3'
 * '<S21>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 4'
 * '<S22>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller x/MRFT Subsystem'
 * '<S23>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller x/MRFT Subsystem/MaxPeak'
 * '<S24>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller x/MRFT Subsystem/MinPeak'
 * '<S25>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller z/MRFT Subsystem'
 * '<S26>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller z/MRFT Subsystem/MaxPeak'
 * '<S27>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /controller z/MRFT Subsystem/MinPeak'
 * '<S28>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 1/reference correction'
 * '<S29>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 2/reference correction'
 * '<S30>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 3/reference correction'
 * '<S31>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/Outer loop Control /convert to horizon 4/reference correction'
 * '<S32>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/R_I_B constructor'
 * '<S33>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/feedback linearization C'
 * '<S34>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/horizon to inertial'
 * '<S35>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/feedback linearization C/Body frame attitude error1'
 * '<S36>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/feedback linearization C/Radians to Degrees'
 * '<S37>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/feedback linearization C/Radians to Degrees1'
 * '<S38>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/feedback linearization C/attitude mapper '
 * '<S39>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Controller/feedback linearization/horizon to inertial/reference correction'
 * '<S40>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Quadcopter Dynamics1/Aerodynamic effect '
 * '<S41>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Quadcopter Dynamics1/Nonlinear Dynamics'
 * '<S42>'  : 'uavNavigationRL2D_2022a_icra/uavController2/Quadcopter Dynamics1/Nonlinear Dynamics/kinematic differential equation'
 */
#endif                          /* RTW_HEADER_uavNavigationRL2D_2022a_icra_h_ */
