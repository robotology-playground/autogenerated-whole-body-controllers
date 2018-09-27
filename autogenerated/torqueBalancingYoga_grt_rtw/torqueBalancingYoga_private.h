/*
 * torqueBalancingYoga_private.h
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueBalancingYoga".
 *
 * Model version              : 1.3292
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C++ source code generated on : Thu Sep 27 14:38:56 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_torqueBalancingYoga_private_h_
#define RTW_HEADER_torqueBalancingYoga_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "torqueBalancingYoga.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm) (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm) (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val) ((rtm)->Timing.t = (val))
#endif

extern void rt_mldivided4x4(const real_T u0[16], const real_T u1[16], real_T y[16]);
extern void fromImuToHomogeousTransformFCN(const real_T rtu_imu_H_link[16],
                                           const real_T rtu_imu_H_link_0[16],
                                           const real_T rtu_link_H_base[16],
                                           const real_T rtu_inertial_0[12],
                                           const real_T rtu_inertial[12],
                                           const real_T rtu_neck_pos[3],
                                           B_fromImuToHomogeousTransform_T* localB,
                                           P_torqueBalancingYoga_T* torqueBalancingYoga_P);
extern void torqueBalan_MATLABFunction_Init(DW_MATLABFunction_torqueBalan_T* localDW);
extern void torqueBalancingY_MATLABFunction(const real_T rtu_s[16],
                                            B_MATLABFunction_torqueBalanc_T* localB,
                                            DW_MATLABFunction_torqueBalan_T* localDW);
extern void torqueBal_MATLABFunction_m_Init(DW_MATLABFunction_torqueBal_k_T* localDW);
extern void torqueBalancin_MATLABFunction_l(const real_T rtu_s[12],
                                            B_MATLABFunction_torqueBala_f_T* localB,
                                            DW_MATLABFunction_torqueBal_k_T* localDW);
extern void torqueBa_MATLABFunction_mb_Init(DW_MATLABFunction_torqueBal_l_T* localDW);
extern void torqueBalancin_MATLABFunction_k(const real_T rtu_s[23],
                                            B_MATLABFunction_torqueBala_j_T* localB,
                                            DW_MATLABFunction_torqueBal_l_T* localDW);
extern void torqueBal_MATLABFunction_l_Init(DW_MATLABFunction_torqueBal_e_T* localDW);
extern void torqueBalancin_MATLABFunction_a(const real_T rtu_s[3],
                                            B_MATLABFunction_torqueBala_p_T* localB,
                                            DW_MATLABFunction_torqueBal_e_T* localDW);

#endif /* RTW_HEADER_torqueBalancingYoga_private_h_ */
