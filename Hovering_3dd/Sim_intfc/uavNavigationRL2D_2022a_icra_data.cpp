/*
 * uavNavigationRL2D_2022a_icra_data.cpp
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

/* Block parameters (default storage) */
P_uavNavigationRL2D_2022a_icr_T uavNavigationRL2D_2022a_icra_P{
  /* Expression: -0.73
   * Referenced by: '<S25>/Beta'
   */
  -0.73,

  /* Expression: 3
   * Referenced by: '<S25>/h'
   */
  3.0,

  /* Expression: 0
   * Referenced by: '<S17>/Switch'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S15>/Constant'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S15>/Constant1'
   */
  0.0,

  /* Expression: -0.73
   * Referenced by: '<S22>/Beta'
   */
  -0.73,

  /* Expression: 8
   * Referenced by: '<S22>/h'
   */
  8.0,

  /* Expression: 1
   * Referenced by: '<S17>/Constant'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S17>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Integrator7'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Integrator6'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Integrator5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay6'
   */
  0.0,

  /* Expression: 180/pi
   * Referenced by: '<S6>/Gain'
   */
  57.295779513082323,

  /* Expression: -1
   * Referenced by: '<S2>/Gain'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay11'
   */
  0.0,

  /* Expression: 180/pi
   * Referenced by: '<S7>/Gain'
   */
  57.295779513082323,

  /* Computed Parameter: TransferFcn7_A
   * Referenced by: '<S5>/Transfer Fcn7'
   */
  -2.0158943593945331,

  /* Computed Parameter: TransferFcn7_C
   * Referenced by: '<S5>/Transfer Fcn7'
   */
  2.0158943593945331,

  /* Computed Parameter: TransferFcn8_A
   * Referenced by: '<S5>/Transfer Fcn8'
   */
  -2.0158943593945331,

  /* Computed Parameter: TransferFcn8_C
   * Referenced by: '<S5>/Transfer Fcn8'
   */
  2.0158943593945331,

  /* Expression: eye(4)
   * Referenced by: '<S5>/Constant2'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  /* Expression: [ 0.0489   -0.4581    0.5566    5.0787;
     0.0489   -0.4581   -0.5566   -5.0787;
     0.0489    0.4581    0.5566   -5.0787;
     0.0489    0.4581   -0.5566    5.0787]
   * Referenced by: '<S5>/Constant4'
   */
  { 0.0489, 0.0489, 0.0489, 0.0489, -0.4581, -0.4581, 0.4581, 0.4581, 0.5566,
    -0.5566, 0.5566, -0.5566, 5.0787, -5.0787, -5.0787, 5.0787 },

  /* Expression: 0
   * Referenced by: '<S3>/Integrator Limited'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Integrator Limited1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Integrator Limited2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Integrator Limited3'
   */
  0.0,

  /* Expression: [Inf Inf Inf Inf]
   * Referenced by: '<S5>/Maximum  Command Authority1'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [-Inf -Inf -Inf -Inf]
   * Referenced by: '<S5>/Maximum  Command Authority1'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: 1
   * Referenced by: '<S5>/Gain13'
   */
  1.0,

  /* Expression: [1 0 0 0]
   * Referenced by: '<S41>/Memory'
   */
  { 1.0, 0.0, 0.0, 0.0 },

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay14'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay14'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/upper_limit_x'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<Root>/upper_limit_y'
   */
  1.0,

  /* Expression: 5
   * Referenced by: '<Root>/upper_limit_z'
   */
  5.0,

  /* Expression: -1
   * Referenced by: '<Root>/lower_limit_x'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<Root>/lower_limit_y'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<Root>/lower_limit_z'
   */
  0.0,

  /* Expression: [ 0.0489   -0.4581    0.5566    5.0787;
     0.0489   -0.4581   -0.5566   -5.0787;
     0.0489    0.4581    0.5566   -5.0787;
     0.0489    0.4581   -0.5566    5.0787]
   * Referenced by: '<S4>/Constant5'
   */
  { 0.0489, 0.0489, 0.0489, 0.0489, -0.4581, -0.4581, 0.4581, 0.4581, 0.5566,
    -0.5566, 0.5566, -0.5566, 5.0787, -5.0787, -5.0787, 5.0787 },

  /* Expression: -1
   * Referenced by: '<S10>/Gain'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<Root>/True1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S2>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S23>/Constant2'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S22>/Gain5'
   */
  -1.0,

  /* Computed Parameter: TransferFcn1_A
   * Referenced by: '<S22>/Transfer Fcn1'
   */
  -20.0,

  /* Computed Parameter: TransferFcn1_C
   * Referenced by: '<S22>/Transfer Fcn1'
   */
  -400.0,

  /* Computed Parameter: TransferFcn1_D
   * Referenced by: '<S22>/Transfer Fcn1'
   */
  20.0,

  /* Expression: 0
   * Referenced by: '<S22>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S23>/Memory'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S22>/Gain4'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S24>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S24>/Memory1'
   */
  0.0,

  /* Expression: 5
   * Referenced by: '<S15>/Switch1'
   */
  5.0,

  /* Expression: 0
   * Referenced by: '<S15>/Switch'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/True2'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/False1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S26>/Constant2'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S25>/Gain5'
   */
  -1.0,

  /* Computed Parameter: TransferFcn1_A_d
   * Referenced by: '<S25>/Transfer Fcn1'
   */
  -20.0,

  /* Computed Parameter: TransferFcn1_C_b
   * Referenced by: '<S25>/Transfer Fcn1'
   */
  -400.0,

  /* Computed Parameter: TransferFcn1_D_d
   * Referenced by: '<S25>/Transfer Fcn1'
   */
  20.0,

  /* Expression: 0
   * Referenced by: '<S25>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S26>/Memory'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S25>/Gain4'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S27>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S27>/Memory1'
   */
  0.0,

  /* Expression: 5
   * Referenced by: '<S17>/Switch1'
   */
  5.0,

  /* Expression: 0
   * Referenced by: '<Root>/False2'
   */
  0.0,

  /* Expression: 0.5
   * Referenced by: '<S9>/Switch'
   */
  0.5,

  /* Expression: 1
   * Referenced by: '<Root>/True3'
   */
  1.0,

  /* Expression: 15
   * Referenced by: '<S2>/max_tilt_ref_limit_deg'
   */
  15.0,

  /* Expression: 9.81
   * Referenced by: '<S33>/virtual anti-gravity (just to make attitude defined)'
   */
  9.81,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay15'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S13>/Gain3'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S5>/Transport Delay16'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S12>/Gain3'
   */
  -1.0,

  /* Expression: [Inf Inf Inf Inf]
   * Referenced by: '<S4>/Maximum  Command Authority'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [-inf -Inf -Inf -Inf]
   * Referenced by: '<S4>/Maximum  Command Authority'
   */
  { 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [max_u max_u max_u max_u]
   * Referenced by: '<S4>/Maximum  Command Authority 1'
   */
  { 1.0, 1.0, 1.0, 1.0 },

  /* Expression: [min_u min_u min_u min_u]
   * Referenced by: '<S4>/Maximum  Command Authority 1'
   */
  { -1.0, -1.0, -1.0, -1.0 },

  /* Expression: 0
   * Referenced by: '<S3>/Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Transport Delay1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Transport Delay2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/Transport Delay3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S40>/Integrator5'
   */
  0.0,

  /* Expression: [0, 0, 1/T_aero]
   * Referenced by: '<S40>/Gain7'
   */
  { 0.0, 0.0, 1.7749773935772943 },

  /* Expression: 0
   * Referenced by: '<S41>/Constant2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S41>/Constant1'
   */
  0.0,

  /* Expression: -9.81
   * Referenced by: '<S5>/Constant3'
   */
  -9.81,

  /* Computed Parameter: Memory2_InitialCondition
   * Referenced by: '<S22>/Memory2'
   */
  false,

  /* Computed Parameter: Memory2_InitialCondition_g
   * Referenced by: '<S25>/Memory2'
   */
  false
};
