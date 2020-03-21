/*
 * torqueControlBalancing_private.h
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueControlBalancing".
 *
 * Model version              : 1.3504
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C++ source code generated on : Sat Mar 21 21:12:50 2020
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_torqueControlBalancing_private_h_
#define RTW_HEADER_torqueControlBalancing_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "torqueControlBalancing.h"

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
extern void rt_invd4x4(const real_T u[16], real_T y[16]);
extern void torqueContr_MATLABFunction_Init(DW_MATLABFunction_torqueContr_T* localDW);
extern void torqueControlBal_MATLABFunction(const real_T rtu_s[23],
                                            B_MATLABFunction_torqueContro_T* localB,
                                            DW_MATLABFunction_torqueContr_T* localDW);
extern void torqueCo_GetBaseRotationFromIMU(const real_T rtu_imu_H_fixedLink[16],
                                            const real_T rtu_imu_H_fixedLink_0[16],
                                            const real_T rtu_fixedLink_H_base[16],
                                            const real_T rtu_rpyFromIMU_0[3],
                                            const real_T rtu_rpyFromIMU[3],
                                            const real_T rtu_neck_pos[3],
                                            B_GetBaseRotationFromIMU_torq_T* localB,
                                            P_torqueControlBalancing_T* torqueControlBalancing_P);
extern void torqueCon_MATLABFunction_i_Init(DW_MATLABFunction_torqueCon_m_T* localDW);
extern void torqueControlB_MATLABFunction_i(const real_T rtu_s[16],
                                            B_MATLABFunction_torqueCont_m_T* localB,
                                            DW_MATLABFunction_torqueCon_m_T* localDW);
extern void torqueCon_MATLABFunction_g_Init(DW_MATLABFunction_torqueCon_p_T* localDW);
extern void torqueControlB_MATLABFunction_o(const real_T rtu_s[12],
                                            B_MATLABFunction_torqueCont_p_T* localB,
                                            DW_MATLABFunction_torqueCon_p_T* localDW);
extern void torqueControlBa_MATLABFunction1(const real_T rtu_d[3],
                                            B_MATLABFunction1_torqueContr_T* localB);

#endif /* RTW_HEADER_torqueControlBalancing_private_h_ */
