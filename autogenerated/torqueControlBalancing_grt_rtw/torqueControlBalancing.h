/*
 * torqueControlBalancing.h
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

#ifndef RTW_HEADER_torqueControlBalancing_h_
#define RTW_HEADER_torqueControlBalancing_h_
#include <cmath>
#include <math.h>
#include <string.h>
#ifndef torqueControlBalancing_COMMON_INCLUDES_
#define torqueControlBalancing_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include <BlockFactory/Core/Block.h>
#include <BlockFactory/Core/FactorySingleton.h>
#include <BlockFactory/Core/Log.h>
#include <BlockFactory/Core/Parameter.h>
#include <BlockFactory/Core/Parameters.h>
#include <BlockFactory/SimulinkCoder/CoderBlockInformation.h>
#include <cstdio>
#endif /* torqueControlBalancing_COMMON_INCLUDES_ */

#include "torqueControlBalancing_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_assert.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm) ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm) (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm) ((rtm)->Timing.t)
#endif

/* Block signals for system '<S59>/MATLAB Function' */
typedef struct
{
    real_T s0[23]; /* '<S59>/MATLAB Function' */
} B_MATLABFunction_torqueContro_T;

/* Block states (default storage) for system '<S59>/MATLAB Function' */
typedef struct
{
    real_T state[23]; /* '<S59>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S59>/MATLAB Function' */
} DW_MATLABFunction_torqueContr_T;

/* Block signals for system '<S70>/Get Base Rotation From IMU' */
typedef struct
{
    real_T w_H_b[16]; /* '<S70>/Get Base Rotation From IMU' */
} B_GetBaseRotationFromIMU_torq_T;

/* Block signals for system '<S77>/MATLAB Function' */
typedef struct
{
    real_T s0[16]; /* '<S77>/MATLAB Function' */
} B_MATLABFunction_torqueCont_m_T;

/* Block states (default storage) for system '<S77>/MATLAB Function' */
typedef struct
{
    real_T state[16]; /* '<S77>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S77>/MATLAB Function' */
} DW_MATLABFunction_torqueCon_m_T;

/* Block signals for system '<S78>/MATLAB Function' */
typedef struct
{
    real_T s0[12]; /* '<S78>/MATLAB Function' */
} B_MATLABFunction_torqueCont_p_T;

/* Block states (default storage) for system '<S78>/MATLAB Function' */
typedef struct
{
    real_T state[12]; /* '<S78>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S78>/MATLAB Function' */
} DW_MATLABFunction_torqueCon_p_T;

/* Block signals for system '<S105>/MATLAB Function1' */
typedef struct
{
    real_T D[9]; /* '<S105>/MATLAB Function1' */
} B_MATLABFunction1_torqueContr_T;

/* Block signals (default storage) */
typedef struct
{
    real_T M_with_inertia[841]; /* '<S45>/Add motor reflected inertias' */
    real_T b_A[841];
    real_T dv0[841];
    real_T pinvLambda[667];
    real_T NullLambda[529];
    real_T NullLambda_Mbar[529];
    real_T D[529]; /* '<S105>/MATLAB Function' */
    real_T dv1[529];
    real_T dv2[529];
    real_T rtb_M_with_inertia_m[529];
    real_T rtb_D_c[529];
    real_T rtb_D_k[529];
    real_T dv3[529];
    real_T dv4[529];
    real_T invTGamma[529];
    real_T invTGamma_t[529];
    real_T T[529];
    real_T b_A_c[529];
    real_T V[529];
    real_T U[529];
    real_T b_A_b[529];
    real_T Vf[529];
    real_T reflectedInertia_p[529];
    real_T invTGamma_c[529];
    real_T invTGamma_t_f[529];
    real_T T_g[529];
    real_T b_A_g[529];
    real_T SFunction[23]; /* '<S3>/S-Function' */
    real_T SFunction_d[23]; /* '<S4>/S-Function' */
    real_T SFunction_k[16]; /* '<S74>/S-Function' */
    real_T SFunction_a[16]; /* '<S72>/S-Function' */
    real_T SFunction_n[16]; /* '<S75>/S-Function' */
    real_T IMU_meas[12]; /* '<Root>/IMU_meas' */
    real_T NeckPosition[3]; /* '<S70>/Neck Position' */
    real_T Switch6[16]; /* '<S70>/Switch6' */
    real_T SFunction_j[16]; /* '<S103>/S-Function' */
    real_T SFunction_f[16]; /* '<S84>/S-Function' */
    real_T SFunction_i[16]; /* '<S73>/S-Function' */
    real_T SFunction_e[16]; /* '<S85>/S-Function' */
    real_T NeckPosition_n[3]; /* '<S71>/Neck Position' */
    real_T Switch6_e[16]; /* '<S71>/Switch6' */
    real_T SFunction_c[16]; /* '<S104>/S-Function' */
    real_T wrench_rightFoot[6]; /* '<Root>/wrench_rightFoot' */
    real_T wrench_leftFoot[6]; /* '<Root>/wrench_leftFoot' */
    real_T MinimumJerkTrajectoryGenerator1[23]; /* '<S107>/MinimumJerkTrajectoryGenerator1' */
    real_T MinimumJerkTrajectoryGenerato_o[23]; /* '<S107>/MinimumJerkTrajectoryGenerator1' */
    real_T MinimumJerkTrajectoryGenerato_p[23]; /* '<S107>/MinimumJerkTrajectoryGenerator1' */
    real_T Switch5[23]; /* '<S107>/Switch5' */
    real_T SFunction_e0[174]; /* '<S68>/S-Function' */
    real_T SFunction_da[174]; /* '<S69>/S-Function' */
    real_T SFunction_nr[23]; /* '<S5>/S-Function' */
    real_T SFunction_b[841]; /* '<S48>/S-Function' */
    real_T Gain[23]; /* '<S47>/Gain' */
    real_T SFunction_ft[29]; /* '<S47>/S-Function' */
    real_T SFunction_ci[6]; /* '<S46>/S-Function' */
    real_T SFunction_br[16]; /* '<S40>/S-Function' */
    real_T SFunction_ai[16]; /* '<S41>/S-Function' */
    real_T Switch[16]; /* '<S21>/Switch' */
    real_T SFunction_e0o[174]; /* '<S37>/S-Function' */
    real_T SFunction_cf[174]; /* '<S38>/S-Function' */
    real_T Add[23]; /* '<S21>/Add' */
    real_T SFunction_o[6]; /* '<S34>/S-Function' */
    real_T SFunction_oe[16]; /* '<S56>/S-Function' */
    real_T SFunction_jt[16]; /* '<S57>/S-Function' */
    real_T SFunction_c4[174]; /* '<S54>/S-Function' */
    real_T SFunction_b3[174]; /* '<S55>/S-Function' */
    real_T SFunction_ej[6]; /* '<S51>/S-Function' */
    real_T SFunction_h[6]; /* '<S52>/S-Function' */
    real_T SFunction_dk[16]; /* '<S50>/S-Function' */
    real_T SFunction_n5[174]; /* '<S53>/S-Function' */
    real_T MinimumJerkTrajectoryGenerator2[3]; /* '<S106>/MinimumJerkTrajectoryGenerator2' */
    real_T MinimumJerkTrajectoryGenerato_a[3]; /* '<S106>/MinimumJerkTrajectoryGenerator2' */
    real_T MinimumJerkTrajectoryGenerat_pp[3]; /* '<S106>/MinimumJerkTrajectoryGenerator2' */
    real_T TmpSignalConversionAtMinimumJer[29];
    real_T MinimumJerkTrajectoryGenerator[29]; /* '<S64>/MinimumJerkTrajectoryGenerator' */
    real_T Saturation[23]; /* '<S17>/Saturation' */
    real_T YarpClock; /* '<S8>/Yarp Clock' */
    real_T noSpikes; /* '<S112>/MATLAB Function' */
    real_T SFunction_o1[23]; /* '<S113>/S-Function' */
    real_T SFunction_o2[23]; /* '<S113>/S-Function' */
    real_T inRange; /* '<S111>/MATLAB Function' */
    real_T w_H_b[16]; /* '<S64>/STATE MACHINE' */
    real_T pos_CoM_des[3]; /* '<S64>/STATE MACHINE' */
    real_T jointPos_des[23]; /* '<S64>/STATE MACHINE' */
    real_T smoothingTimeJoints; /* '<S64>/STATE MACHINE' */
    real_T smoothingTimeCoM; /* '<S64>/STATE MACHINE' */
    real_T nu_b[6]; /* '<S62>/Compute Base Velocity' */
    real_T HessianMatrixOneFoot[36]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T gradientOneFoot[6]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T ConstraintsMatrixOneFoot[114]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T bVectorConstraintsOneFoot[19]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T HessianMatrixTwoFeet[144]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T gradientTwoFeet[12]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T ConstraintsMatrixTwoFeet[456]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T bVectorConstraintsTwoFeet[38]; /* '<S15>/Momentum Based Balancing Controller ' */
    real_T reflectedInertia[529]; /* '<S22>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}' */
    real_T baseVel_equivalent[6]; /* '<S21>/Get Equivalent Base Velocity' */
    real_T QPTwoFeet_o1[12]; /* '<S27>/QP Two Feet' */
    real_T QPTwoFeet_o2; /* '<S27>/QP Two Feet' */
    real_T f_star[12]; /* '<S27>/Process QP output' */
    real_T QPOneFoot_o1[6]; /* '<S26>/QP One Foot' */
    real_T QPOneFoot_o2; /* '<S26>/QP One Foot' */
    real_T f_star_m[12]; /* '<S26>/Process QP output' */
    real_T SFunction_bm[23]; /* '<S14>/S-Function' */
    boolean_T HiddenBuf_InsertedFor_QPOneFoot; /* '<S24>/One Foot Two Feet QP Selector' */
    boolean_T NOT; /* '<S24>/NOT' */
    boolean_T HiddenBuf_InsertedFor_QPTwoFeet; /* '<S24>/NOT' */
    boolean_T onOneFoot; /* '<S24>/One Foot Two Feet QP Selector' */
    B_MATLABFunction1_torqueContr_T sf_MATLABFunction2; /* '<S105>/MATLAB Function2' */
    B_MATLABFunction1_torqueContr_T sf_MATLABFunction1; /* '<S105>/MATLAB Function1' */
    B_MATLABFunction_torqueContro_T sf_MATLABFunction_d; /* '<S96>/MATLAB Function' */
    B_MATLABFunction_torqueCont_p_T sf_MATLABFunction_c; /* '<S88>/MATLAB Function' */
    B_MATLABFunction_torqueCont_m_T sf_MATLABFunction_j; /* '<S87>/MATLAB Function' */
    B_GetBaseRotationFromIMU_torq_T
        sf_GetBaseRotationFromIMU_c; /* '<S71>/Get Base Rotation From IMU' */
    B_MATLABFunction_torqueCont_p_T sf_MATLABFunction_o; /* '<S78>/MATLAB Function' */
    B_MATLABFunction_torqueCont_m_T sf_MATLABFunction_i; /* '<S77>/MATLAB Function' */
    B_GetBaseRotationFromIMU_torq_T
        sf_GetBaseRotationFromIMU; /* '<S70>/Get Base Rotation From IMU' */
    B_MATLABFunction_torqueContro_T sf_MATLABFunction; /* '<S59>/MATLAB Function' */
} B_torqueControlBalancing_T;

/* Block states (default storage) for system '<Root>' */
typedef struct
{
    real_T UnitDelay_DSTATE; /* '<S112>/Unit Delay' */
    real_T UnitDelay_DSTATE_f; /* '<S111>/Unit Delay' */
    real_T u_previous[23]; /* '<S112>/MATLAB Function' */
    real_T state[3]; /* '<S95>/MATLAB Function' */
    real_T currentState; /* '<S64>/STATE MACHINE' */
    real_T t_switch; /* '<S64>/STATE MACHINE' */
    real_T w_H_fixedLink[16]; /* '<S64>/STATE MACHINE' */
    real_T yogaMovesetCounter; /* '<S64>/STATE MACHINE' */
    real_T uPrev[23]; /* '<S17>/Saturate Torque Derivative' */
    real_T inv_DWORK4[16]; /* '<S35>/inv  ' */
    real_T inv1_DWORK4[16]; /* '<S35>/inv  1' */
    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK; /* '<S3>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_c; /* '<S4>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_p; /* '<S74>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_b; /* '<S72>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_m; /* '<S75>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } IMU_meas_PWORK; /* '<Root>/IMU_meas' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK; /* '<S70>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_h; /* '<S103>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_g; /* '<S84>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_l; /* '<S73>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_a; /* '<S85>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK_j; /* '<S71>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_k; /* '<S104>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } wrench_rightFoot_PWORK; /* '<Root>/wrench_rightFoot' */

    struct
    {
        void* blockPWork[2];
    } wrench_leftFoot_PWORK; /* '<Root>/wrench_leftFoot' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator1; /* '<S107>/MinimumJerkTrajectoryGenerator1' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_j; /* '<S68>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_jg; /* '<S69>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_bs; /* '<S5>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_n; /* '<S48>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_hk; /* '<S47>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_e; /* '<S46>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d; /* '<S40>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_er; /* '<S41>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_dc; /* '<S37>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_am; /* '<S38>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_c4; /* '<S34>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_jj; /* '<S56>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_pk; /* '<S57>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_ne; /* '<S54>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_bp; /* '<S55>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d4; /* '<S51>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d1; /* '<S52>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_hl; /* '<S50>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_ay; /* '<S53>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator2; /* '<S106>/MinimumJerkTrajectoryGenerator2' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator_; /* '<S64>/MinimumJerkTrajectoryGenerator' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_le; /* '<S7>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } YarpClock_PWORK; /* '<S8>/Yarp Clock' */

    struct
    {
        void* blockPWork[2];
    } RealTimeSynchronizer_PWORK; /* '<S117>/Real Time Synchronizer' */

    struct
    {
        void* blockPWork[2];
    } SimulatorSynchronizer_PWORK; /* '<S116>/Simulator Synchronizer' */

    void* Assertion_slioAccessor; /* '<S112>/Assertion' */
    void* Assertion_slioAccessor_k; /* '<S111>/Assertion' */
    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_g2; /* '<S113>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } QPTwoFeet_PWORK; /* '<S27>/QP Two Feet' */

    struct
    {
        void* blockPWork[2];
    } QPOneFoot_PWORK; /* '<S26>/QP One Foot' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_f; /* '<S14>/S-Function' */

    boolean_T u_previous_not_empty; /* '<S112>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S95>/MATLAB Function' */
    boolean_T currentState_not_empty; /* '<S64>/STATE MACHINE' */
    boolean_T uPrev_not_empty; /* '<S17>/Saturate Torque Derivative' */
    boolean_T
        STOPIFTHEREARESPIKESINTHEENCODE; /* '<S19>/STOP IF THERE ARE SPIKES IN THE ENCODERS' */
    boolean_T STOPIFJOINTSHITTHELIMITS_MODE; /* '<S19>/STOP IF JOINTS HIT THE LIMITS' */
    boolean_T QPTwoFeet_MODE; /* '<S24>/QP Two Feet' */
    boolean_T QPOneFoot_MODE; /* '<S24>/QP One Foot' */
    boolean_T Visualizer_MODE; /* '<S2>/Visualizer' */
    boolean_T DesiredandMeasuredForces_MODE; /* '<S2>/Desired and Measured Forces' */
    DW_MATLABFunction_torqueContr_T sf_MATLABFunction_d; /* '<S96>/MATLAB Function' */
    DW_MATLABFunction_torqueCon_p_T sf_MATLABFunction_c; /* '<S88>/MATLAB Function' */
    DW_MATLABFunction_torqueCon_m_T sf_MATLABFunction_j; /* '<S87>/MATLAB Function' */
    DW_MATLABFunction_torqueCon_p_T sf_MATLABFunction_o; /* '<S78>/MATLAB Function' */
    DW_MATLABFunction_torqueCon_m_T sf_MATLABFunction_i; /* '<S77>/MATLAB Function' */
    DW_MATLABFunction_torqueContr_T sf_MATLABFunction; /* '<S59>/MATLAB Function' */
} DW_torqueControlBalancing_T;

/* Parameters (default storage) */
struct P_torqueControlBalancing_T_
{
    struct_cik6KQn5YVVa3skh1Wn00 StateMachine; /* Variable: StateMachine
                                                * Referenced by: '<S64>/STATE MACHINE'
                                                */
    struct_8fk7vJtSw4UnB4xn6cTd0B Config; /* Variable: Config
                                           * Referenced by:
                                           *   '<S2>/ON_GAZEBO 1'
                                           *   '<S2>/ON_GAZEBO 2'
                                           *   '<S2>/ON_GAZEBO 4'
                                           *   '<S2>/ON_GAZEBO 6'
                                           *   '<S8>/ON_GAZEBO '
                                           *   '<S15>/Momentum Based Balancing Controller '
                                           *   '<S17>/Constant'
                                           *   '<S17>/Constant1'
                                           *   '<S19>/ON_GAZEBO 1'
                                           *   '<S19>/ON_GAZEBO 2'
                                           *   '<S22>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}'
                                           *   '<S22>/ '
                                           *   '<S22>/ 1'
                                           *   '<S22>/Gain'
                                           *   '<S64>/STATE MACHINE'
                                           *   '<S45>/Add motor reflected inertias'
                                           *   '<S45>/Constant'
                                           *   '<S70>/Get Base Rotation From IMU'
                                           *   '<S70>/USE_IMU4EST_BASE'
                                           *   '<S71>/Get Base Rotation From IMU'
                                           *   '<S71>/USE_IMU4EST_BASE'
                                           *   '<S106>/SMOOTH_DES_COM'
                                           *   '<S107>/SMOOTH_DES_COM2'
                                           *   '<S26>/Process QP output'
                                           *   '<S27>/Process QP output'
                                           *   '<S79>/USE_IMU4EST_BASE1'
                                           *   '<S89>/USE_IMU4EST_BASE1'
                                           */
    struct_glrlyr9IieN6zDechOWztB Gain; /* Variable: Gain
                                         * Referenced by:
                                         *   '<S15>/Momentum Based Balancing Controller '
                                         *   '<S64>/STATE MACHINE'
                                         */
    struct_Yt0Pn6XEWtmZM9M9GtOwoD Reg; /* Variable: Reg
                                        * Referenced by:
                                        *   '<S15>/Momentum Based Balancing Controller '
                                        *   '<S21>/Get Equivalent Base Velocity'
                                        *   '<S62>/Compute Base Velocity'
                                        */
    struct_8ovnsQmVgsrxBHMDpqGPN Sat; /* Variable: Sat
                                       * Referenced by:
                                       *   '<S17>/Constant2'
                                       *   '<S17>/Saturation'
                                       *   '<S112>/index1'
                                       */
    real_T ConstraintsMatrix[114]; /* Variable: ConstraintsMatrix
                                    * Referenced by: '<S15>/Constant'
                                    */
    real_T bVectorConstraints[19]; /* Variable: bVectorConstraints
                                    * Referenced by: '<S15>/Constant1'
                                    */
    real_T CompareToConstant_const; /* Mask Parameter: CompareToConstant_const
                                     * Referenced by: '<S80>/Constant'
                                     */
    real_T CompareToConstant_const_o; /* Mask Parameter: CompareToConstant_const_o
                                       * Referenced by: '<S82>/Constant'
                                       */
    real_T CompareToConstant_const_p; /* Mask Parameter: CompareToConstant_const_p
                                       * Referenced by: '<S99>/Constant'
                                       */
    real_T CompareToConstant_const_pf; /* Mask Parameter: CompareToConstant_const_pf
                                        * Referenced by: '<S101>/Constant'
                                        */
    real_T CompareToConstant_const_h; /* Mask Parameter: CompareToConstant_const_h
                                       * Referenced by: '<S90>/Constant'
                                       */
    real_T CompareToConstant_const_k; /* Mask Parameter: CompareToConstant_const_k
                                       * Referenced by: '<S92>/Constant'
                                       */
    real_T CompareToConstant_const_l; /* Mask Parameter: CompareToConstant_const_l
                                       * Referenced by: '<S60>/Constant'
                                       */
    real_T Gain_Gain; /* Expression: pi/180
                       * Referenced by: '<S79>/Gain'
                       */
    real_T Gain_Gain_j; /* Expression: pi/180
                         * Referenced by: '<S89>/Gain'
                         */
    real_T UnitDelay_InitialCondition; /* Expression: 1
                                        * Referenced by: '<S111>/Unit Delay'
                                        */
    real_T index1_Value; /* Expression: 0.01
                          * Referenced by: '<S111>/index1'
                          */
    real_T UnitDelay_InitialCondition_j; /* Expression: 1
                                          * Referenced by: '<S112>/Unit Delay'
                                          */
    real_T Constant7_Value[16]; /* Expression: eye(4)
                                 * Referenced by: '<S63>/Constant7'
                                 */
    real_T Constant_Value[3]; /* Expression: zeros(3,1)
                               * Referenced by: '<S79>/Constant'
                               */
    real_T Constant_Value_e[3]; /* Expression: zeros(3,1)
                                 * Referenced by: '<S89>/Constant'
                                 */
    real_T Constant_Value_d[6]; /* Expression: zeros(6,1)
                                 * Referenced by: '<S47>/Constant'
                                 */
    real_T Gain_Gain_o; /* Expression: 0
                         * Referenced by: '<S47>/Gain'
                         */
    real_T Constant7_Value_o[16]; /* Expression: eye(4)
                                   * Referenced by: '<S35>/Constant7'
                                   */
    real_T Switch_Threshold; /* Expression: 0.1
                              * Referenced by: '<S21>/Switch'
                              */
    real_T Constant_Value_i[6]; /* Expression: zeros(6,1)
                                 * Referenced by: '<S106>/Constant'
                                 */
};

/* Real-time Model Data Structure */
struct tag_RTM_torqueControlBalancin_T
{
    const char_T* errorStatus;
    RTWSolverInfo solverInfo;

    /*
     * Timing:
     * The following substructure contains information regarding
     * the timing information for the model.
     */
    struct
    {
        uint32_T clockTick0;
        uint32_T clockTickH0;
        time_T stepSize0;
        uint32_T clockTick1;
        uint32_T clockTickH1;
        SimTimeStep simTimeStep;
        time_T* t;
        time_T tArray[2];
    } Timing;
};

/* Class declaration for model torqueControlBalancing */
class torqueControlBalancingModelClass
{
    /* public data and function members */
public:
    /* model initialize function */
    void initialize();

    /* model step function */
    void step();

    /* model terminate function */
    void terminate();

    /* Constructor */
    torqueControlBalancingModelClass();

    /* Destructor */
    ~torqueControlBalancingModelClass();

    /* Real-Time Model get method */
    RT_MODEL_torqueControlBalanci_T* getRTM();

    /* private data and function members */
private:
    /* Tunable parameters */
    P_torqueControlBalancing_T torqueControlBalancing_P;

    /* Block signals */
    B_torqueControlBalancing_T torqueControlBalancing_B;

    /* Block states */
    DW_torqueControlBalancing_T torqueControlBalancing_DW;

    /* Real-Time Model */
    RT_MODEL_torqueControlBalanci_T torqueControlBalancing_M;

    /* private member function(s) for subsystem '<Root>'*/
    void torqueControlBalancing_xswap(real_T x[529], int32_T ix0, int32_T iy0);
    void torqueControlBalancing_xgetrf_e(real_T A[529], int32_T ipiv[23], int32_T* info);
    void torqueControlBalancing_xtrsm(const real_T A[529], real_T B[529]);
    void torqueControlBalancing_xtrsm_i(const real_T A[529], real_T B[529]);
    void torqueControlBalancing_mrdivide(const real_T A[529], const real_T B[529], real_T y[529]);
    void t_computeMotorsReflectedInertia(const real_T Gamma[529],
                                         const real_T T[529],
                                         const real_T I_m[529],
                                         real_T reflectedInertia[529]);
    void torqueControlBalancing_xgetrf_a(real_T A[16], int32_T ipiv[4], int32_T* info);
    void torqueControlBalanci_mrdivide_h(real_T A[16], const real_T B[16]);
    void torqueControlBalanci_mldivide_o(const real_T A[16], real_T B[4]);
    real_T torqueControlBalancing_norm(const real_T x[6]);
    void torqueControlBalan_pinvDamped_i(const real_T A[72], real_T regDamp, real_T pinvDampA[72]);
    void torqu_addMotorsReflectedInertia(const real_T M[841],
                                         const real_T Gamma[529],
                                         const real_T T[529],
                                         const real_T I_m[529],
                                         real_T M_with_inertia[841]);
    void torqueControlBalanci_mldivide_l(const real_T A[36], real_T B[36]);
    void torqueControlBalancing_blkdiag(const real_T varargin_1[9],
                                        const real_T varargin_2[9],
                                        real_T y[36]);
    real_T torqueControlBalancing_xnrm2(int32_T n, const real_T x[72], int32_T ix0);
    real_T torqueControlBalancing_xnrm2_h(int32_T n, const real_T x[6], int32_T ix0);
    void torqueControlBalancing_xaxpy_cf(int32_T n,
                                         real_T a,
                                         const real_T x[12],
                                         int32_T ix0,
                                         real_T y[72],
                                         int32_T iy0);
    void torqueControlBalancing_xaxpy_c(int32_T n,
                                        real_T a,
                                        const real_T x[72],
                                        int32_T ix0,
                                        real_T y[12],
                                        int32_T iy0);
    real_T torqueControlBalancing_xdotc(int32_T n,
                                        const real_T x[72],
                                        int32_T ix0,
                                        const real_T y[72],
                                        int32_T iy0);
    void torqueControlBalancing_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[72], int32_T iy0);
    real_T torqueControlBalancing_xdotc_i(int32_T n,
                                          const real_T x[36],
                                          int32_T ix0,
                                          const real_T y[36],
                                          int32_T iy0);
    void
    torqueControlBalancin_xaxpy_cfu(int32_T n, real_T a, int32_T ix0, real_T y[36], int32_T iy0);
    void torqueControlBalancing_xscal(real_T a, real_T x[72], int32_T ix0);
    void torqueControlBalancing_xscal_m(real_T a, real_T x[36], int32_T ix0);
    void torqueControlBalancing_xswap_h(real_T x[36], int32_T ix0, int32_T iy0);
    void torqueControlBalancing_xswap_hc(real_T x[72], int32_T ix0, int32_T iy0);
    void torqueControlBalancing_xrotg(real_T* a, real_T* b, real_T* c, real_T* s);
    void torqueControlBalancing_xrot(real_T x[36], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void torqueControlBalancing_xrot_f(real_T x[72], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void torqueControlBalancing_svd(const real_T A[72], real_T U[72], real_T s[6], real_T V[36]);
    void torqueControlBalancing_pinv(const real_T A[72], real_T tol, real_T X[72]);
    void torqueControlBalancing_eye(real_T I[144]);
    void torqueControlBalancin_xswap_hco(real_T x[841], int32_T ix0, int32_T iy0);
    void torqueControlBalancing_xgetrf_g(real_T A[841], int32_T ipiv[29], int32_T* info);
    void torqueControlBalancing_xtrsm_j(const real_T A[841], real_T B[348]);
    void torqueControlBalancing_xtrsm_jj(const real_T A[841], real_T B[348]);
    void torqueControlBalanci_mrdivide_b(const real_T A[348], const real_T B[841], real_T y[348]);
    void
    torqueControlBalanci_pinvDamped(const real_T A[276], real_T regDamp, real_T pinvDampA[276]);
    void torqueControlBalancing_eye_j(real_T I[529]);
    void torqueControlBalanc_mrdivide_bo(real_T A[138], const real_T B[36]);
    real_T torqueControlBalancing_xnrm2_ht(int32_T n, const real_T x[529], int32_T ix0);
    real_T torqueControlBalancin_xnrm2_ht1(int32_T n, const real_T x[23], int32_T ix0);
    void torqueControlBalan_xaxpy_cfudir(int32_T n,
                                         real_T a,
                                         const real_T x[23],
                                         int32_T ix0,
                                         real_T y[529],
                                         int32_T iy0);
    void torqueControlBalanc_xaxpy_cfudi(int32_T n,
                                         real_T a,
                                         const real_T x[529],
                                         int32_T ix0,
                                         real_T y[23],
                                         int32_T iy0);
    real_T torqueControlBalancing_xdotc_i5(int32_T n,
                                           const real_T x[529],
                                           int32_T ix0,
                                           const real_T y[529],
                                           int32_T iy0);
    void
    torqueControlBalanci_xaxpy_cfud(int32_T n, real_T a, int32_T ix0, real_T y[529], int32_T iy0);
    void torqueControlBalancing_xscal_ms(real_T a, real_T x[529], int32_T ix0);
    void torqueControlBalanci_xswap_hcop(real_T x[529], int32_T ix0, int32_T iy0);
    void
    torqueControlBalancing_xrot_fl(real_T x[529], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void
    torqueControlBalancing_svd_g(const real_T A[529], real_T U[529], real_T s[23], real_T V[529]);
    void torqueControlBalancing_xgemm(int32_T k,
                                      const real_T A[529],
                                      const real_T B[529],
                                      real_T C[529]);
    void torqueControlBalancing_pinv_o(const real_T A[529], real_T tol, real_T X[529]);
    void torqueControlBalancing_diag(const real_T v[23], real_T d[529]);
    void torqueControlBalancin_blkdiag_h(const real_T varargin_1[114],
                                         const real_T varargin_2[114],
                                         real_T y[456]);
    void torqueControlBalancing_invNxN(const real_T x[36], real_T y[36]);
    void torqueControlBalancing_invNxN_b(const real_T x[144], real_T y[144]);
};

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
 * '<Root>' : 'torqueControlBalancing'
 * '<S1>'   : 'torqueControlBalancing/Configuration'
 * '<S2>'   : 'torqueControlBalancing/Dump and visualize'
 * '<S3>'   : 'torqueControlBalancing/GetMeasurement'
 * '<S4>'   : 'torqueControlBalancing/GetMeasurement1'
 * '<S5>'   : 'torqueControlBalancing/GetMeasurement2'
 * '<S6>'   : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL'
 * '<S7>'   : 'torqueControlBalancing/SetReferences'
 * '<S8>'   : 'torqueControlBalancing/synchronizer'
 * '<S9>'   : 'torqueControlBalancing/Dump and visualize/Desired and Measured Forces'
 * '<S10>'  : 'torqueControlBalancing/Dump and visualize/Feet Status and Gains'
 * '<S11>'  : 'torqueControlBalancing/Dump and visualize/Visualize eventual QP failures'
 * '<S12>'  : 'torqueControlBalancing/Dump and visualize/Visualizer'
 * '<S13>'  : 'torqueControlBalancing/Dump and visualize/Desired and Measured Forces/Demux Forces nd
 * Moments'
 * '<S14>'  : 'torqueControlBalancing/Dump and visualize/Visualizer/GetMeasurement1'
 * '<S15>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP'
 * '<S16>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and Kinematics'
 * '<S17>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Joint Torque Saturation'
 * '<S18>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References'
 * '<S19>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint limits'
 * '<S20>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques'
 * '<S21>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error'
 * '<S22>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/From
 * tau_QP to Joint Torques (motor reflected inertia)'
 * '<S23>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Momentum
 * Based Balancing Controller '
 * '<S24>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver'
 * '<S25>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/One Foot Two Feet QP Selector'
 * '<S26>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP One Foot'
 * '<S27>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP Two Feet'
 * '<S28>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP One Foot/Analytical Solution QP One Foot (unconstrained)'
 * '<S29>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP One Foot/MatchSignalSizes'
 * '<S30>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP One Foot/Process QP output'
 * '<S31>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP Two Feet/Analytical Solution Two Feet (unconstrained)'
 * '<S32>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP Two Feet/MatchSignalSizes'
 * '<S33>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * Desired Torques/QPSolver/QP Two Feet/Process QP output'
 * '<S34>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/CentroidalMomentum'
 * '<S35>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Compute base to fixed link transform'
 * '<S36>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Get Equivalent Base Velocity'
 * '<S37>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Jacobian'
 * '<S38>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Jacobian1'
 * '<S39>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Select base to world transform'
 * '<S40>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Compute base to fixed link transform/l_sole to root_link relative
 * transform'
 * '<S41>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/Compute
 * angular momentum integral error/Compute base to fixed link transform/r_sole to root_link relative
 * transform'
 * '<S42>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Balancing Controller QP/From
 * tau_QP to Joint Torques (motor reflected inertia)/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}'
 * '<S43>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics'
 * '<S44>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics'
 * '<S45>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics/Add motors reflected inertia'
 * '<S46>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics/CentroidalMomentum'
 * '<S47>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics/GetBiasForces'
 * '<S48>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics/MassMatrix'
 * '<S49>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Dynamics/Add motors reflected inertia/Add motor reflected inertias'
 * '<S50>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/CoM'
 * '<S51>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/DotJNu l_sole '
 * '<S52>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/DotJNu r_sole   '
 * '<S53>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/Jacobian com'
 * '<S54>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/Jacobian l_sole'
 * '<S55>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/Jacobian r_sole'
 * '<S56>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/l_sole'
 * '<S57>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Dynamics and
 * Kinematics/Kinematics/r_sole'
 * '<S58>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Joint Torque Saturation/Saturate
 * Torque Derivative'
 * '<S59>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Joint Torque Saturation/holder '
 * '<S60>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Joint Torque Saturation/holder
 * /Compare To Constant'
 * '<S61>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Joint Torque Saturation/holder
 * /MATLAB Function'
 * '<S62>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute State Velocity'
 * '<S63>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform'
 * '<S64>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine'
 * '<S65>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References'
 * '<S66>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute State Velocity/Compute Base Velocity'
 * '<S67>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute State Velocity/Feet Jacobians'
 * '<S68>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute State Velocity/Feet Jacobians/Jacobian l_sole'
 * '<S69>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute State Velocity/Feet Jacobians/Jacobian r_sole'
 * '<S70>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform '
 * '<S71>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform'
 * '<S72>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/Relative transform l_sole_H_base'
 * '<S73>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/Relative transform r_sole_H_base'
 * '<S74>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /Fixed base to imu
 * transform'
 * '<S75>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /Fixed base to root
 * link transform'
 * '<S76>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /Get Base Rotation
 * From IMU'
 * '<S77>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 1'
 * '<S78>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 2'
 * '<S79>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /neck transform'
 * '<S80>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 1/Compare To
 * Constant'
 * '<S81>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 1/MATLAB
 * Function'
 * '<S82>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 2/Compare To
 * Constant'
 * '<S83>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/LFoot to base link transform /holder 2/MATLAB
 * Function'
 * '<S84>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/Fixed base to imu
 * transform'
 * '<S85>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/Fixed base to root
 * link transform'
 * '<S86>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/Get Base Rotation
 * From IMU'
 * '<S87>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 1'
 * '<S88>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 2'
 * '<S89>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/neck transform'
 * '<S90>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 1/Compare To
 * Constant'
 * '<S91>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 1/MATLAB
 * Function'
 * '<S92>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 2/Compare To
 * Constant'
 * '<S93>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Compute base to fixed link transform/RFoot to base link transform/holder 2/MATLAB
 * Function'
 * '<S94>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/STATE MACHINE'
 * '<S95>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 1'
 * '<S96>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 2'
 * '<S97>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/xCom'
 * '<S98>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/xCom2'
 * '<S99>'  : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 1/Compare To Constant'
 * '<S100>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 1/MATLAB Function'
 * '<S101>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 2/Compare To Constant'
 * '<S102>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/holder 2/MATLAB Function'
 * '<S103>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/xCom/CoM'
 * '<S104>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and References/State
 * Machine/xCom2/CoM'
 * '<S105>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Reshape Gains Matrices'
 * '<S106>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Smooth reference CoM'
 * '<S107>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Smooth reference joint position'
 * '<S108>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Reshape Gains Matrices/MATLAB Function'
 * '<S109>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Reshape Gains Matrices/MATLAB Function1'
 * '<S110>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/Robot State and
 * References/Update Gains and References/Reshape Gains Matrices/MATLAB Function2'
 * '<S111>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint
 * limits/STOP IF JOINTS HIT THE LIMITS'
 * '<S112>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint
 * limits/STOP IF THERE ARE SPIKES IN THE ENCODERS'
 * '<S113>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint
 * limits/STOP IF JOINTS HIT THE LIMITS/GetLimits'
 * '<S114>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint
 * limits/STOP IF JOINTS HIT THE LIMITS/MATLAB Function'
 * '<S115>' : 'torqueControlBalancing/MOMENTUM BASED TORQUE CONTROL/emergency stop: joint
 * limits/STOP IF THERE ARE SPIKES IN THE ENCODERS/MATLAB Function'
 * '<S116>' : 'torqueControlBalancing/synchronizer/GAZEBO_SYNCHRONIZER'
 * '<S117>' : 'torqueControlBalancing/synchronizer/REAL_TIME_SYNC'
 */
#endif /* RTW_HEADER_torqueControlBalancing_h_ */
