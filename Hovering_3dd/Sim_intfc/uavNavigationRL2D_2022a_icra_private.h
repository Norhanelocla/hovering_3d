/*
 * uavNavigationRL2D_2022a_icra_private.h
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

#ifndef RTW_HEADER_uavNavigationRL2D_2022a_icra_private_h_
#define RTW_HEADER_uavNavigationRL2D_2022a_icra_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "uavNavigationRL2D_2022a_icra.h"
#include "uavNavigationRL2D_2022a_icra_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

#if !defined(rt_VALIDATE_MEMORY)
#define rt_VALIDATE_MEMORY(S, ptr)     if(!(ptr)) {\
 rtmSetErrorStatus(uavNavigationRL2D_2022a_icra_M, RT_MEMORY_ALLOCATION_ERROR);\
 }
#endif

#if !defined(rt_FREE)
#if !defined(_WIN32)
#define rt_FREE(ptr)                   if((ptr) != (nullptr)) {\
 free((ptr));\
 (ptr) = (nullptr);\
 }
#else

/* Visual and other windows compilers declare free without const */
#define rt_FREE(ptr)                   if((ptr) != (nullptr)) {\
 free((void *)(ptr));\
 (ptr) = (nullptr);\
 }
#endif
#endif

extern void rt_mrdivided4x4_snf(const real_T u0[16], const real_T u1[16], real_T
  y[16]);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
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
  ;
void uavNavigati_referencecorrection(real_T rtu_yaw, real_T rtu_x, real_T rtu_y,
  B_referencecorrection_uavNavi_T *localB);
void uavNav_referencecorrection_Term(void);

/* private model entry point functions */
extern void uavNavigationRL2D_2022a_icra_derivatives
  (RT_MODEL_uavNavigationRL2D_20_T *const uavNavigationRL2D_2022a_icra_M);

#endif                  /* RTW_HEADER_uavNavigationRL2D_2022a_icra_private_h_ */
