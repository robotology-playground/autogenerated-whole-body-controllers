/*
 * torqueBalancingYoga.h
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueBalancingYoga".
 *
 * Model version              : 1.3253
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C++ source code generated on : Wed May 23 14:31:10 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_torqueBalancingYoga_h_
#define RTW_HEADER_torqueBalancingYoga_h_
#include <cmath>
#include <math.h>
#include <string.h>
#ifndef torqueBalancingYoga_COMMON_INCLUDES_
#define torqueBalancingYoga_COMMON_INCLUDES_
#include <cstdio>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "Block.h"
#include "Log.h"
#include "Parameter.h"
#include "Parameters.h"
#include "CoderBlockInformation.h"
#include "ModelPartitioner.h"
#include "GetMeasurement.h"
#include "ForwardKinematics.h"
#include "YarpRead.h"
#include "Jacobian.h"
#include "MinimumJerkTrajectoryGenerator.h"
#include "QpOases.h"
#include "GetLimits.h"
#include "SimulatorSynchronizer.h"
#include "RealTimeSynchronizer.h"
#include "MassMatrix.h"
#include "InverseDynamics.h"
#include "CentroidalMomentum.h"
#include "DotJNu.h"
#include "SetReferences.h"
#endif /* torqueBalancingYoga_COMMON_INCLUDES_ */

#include "torqueBalancingYoga_types.h"

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

/* Block signals for system '<S41>/fromImuToHomogeousTransformFCN' */
typedef struct
{
    real_T w_H_b[16]; /* '<S41>/fromImuToHomogeousTransformFCN' */
} B_fromImuToHomogeousTransform_T;

/* Block signals for system '<S47>/MATLAB Function' */
typedef struct
{
    real_T s0[16]; /* '<S47>/MATLAB Function' */
} B_MATLABFunction_torqueBalanc_T;

/* Block states (auto storage) for system '<S47>/MATLAB Function' */
typedef struct
{
    real_T state[16]; /* '<S47>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S47>/MATLAB Function' */
} DW_MATLABFunction_torqueBalan_T;

/* Block signals for system '<S48>/MATLAB Function' */
typedef struct
{
    real_T s0[12]; /* '<S48>/MATLAB Function' */
} B_MATLABFunction_torqueBala_f_T;

/* Block states (auto storage) for system '<S48>/MATLAB Function' */
typedef struct
{
    real_T state[12]; /* '<S48>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S48>/MATLAB Function' */
} DW_MATLABFunction_torqueBal_k_T;

/* Block signals for system '<S38>/MATLAB Function' */
typedef struct
{
    real_T s0[23]; /* '<S38>/MATLAB Function' */
} B_MATLABFunction_torqueBala_j_T;

/* Block states (auto storage) for system '<S38>/MATLAB Function' */
typedef struct
{
    real_T state[23]; /* '<S38>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S38>/MATLAB Function' */
} DW_MATLABFunction_torqueBal_l_T;

/* Block signals for system '<S39>/MATLAB Function' */
typedef struct
{
    real_T s0[3]; /* '<S39>/MATLAB Function' */
} B_MATLABFunction_torqueBala_p_T;

/* Block states (auto storage) for system '<S39>/MATLAB Function' */
typedef struct
{
    real_T state[3]; /* '<S39>/MATLAB Function' */
    boolean_T state_not_empty; /* '<S39>/MATLAB Function' */
} DW_MATLABFunction_torqueBal_e_T;

/* Block signals (auto storage) */
typedef struct
{
    real_T M_with_inertia[841]; /* '<S15>/Add motor reflected inertias' */
    real_T b_A[841];
    real_T Gamma[529];
    real_T invTGamma[529];
    real_T dv0[529];
    real_T dv1[529];
    real_T dv2[529];
    real_T Gamma_m[529];
    real_T dv3[529];
    real_T dv4[529];
    real_T invTGamma_c[529];
    real_T invTGamma_t[529];
    real_T b_Config[529];
    real_T b_A_k[529];
    real_T b_A_c[529];
    real_T V[529];
    real_T U[529];
    real_T b_A_b[529];
    real_T Vf[529];
    real_T b_A_p[529];
    real_T SFunction[23]; /* '<S151>/S-Function' */
    real_T SFunction_d[23]; /* '<S31>/S-Function' */
    real_T MultiportSwitch1[133]; /* '<S4>/Multiport Switch1' */
    real_T u[529]; /* '<S6>/    5' */
    real_T Constant[114]; /* '<S6>/Constant' */
    real_T Constant1[19]; /* '<S6>/Constant1' */
    real_T MinimumJerkTrajectoryGenerator1[23]; /* '<S33>/Minimum Jerk Trajectory Generator1' */
    real_T Switch5[23]; /* '<S33>/Switch5' */
    real_T SFunction_b[841]; /* '<S18>/S-Function' */
    real_T Constant_h[6]; /* '<S17>/Constant' */
    real_T Gain[23]; /* '<S17>/Gain' */
    real_T SFunction_f[29]; /* '<S17>/S-Function' */
    real_T SFunction_c[6]; /* '<S16>/S-Function' */
    real_T Constant7[16]; /* '<S113>/Constant7' */
    real_T SFunction_bd[16]; /* '<S122>/S-Function' */
    real_T SFunction_br[16]; /* '<S119>/S-Function' */
    real_T Constant_d; /* '<S128>/Constant' */
    real_T SFunction_d2[16]; /* '<S123>/S-Function' */
    real_T inertial[12]; /* '<S109>/inertial' */
    real_T Constant_m; /* '<S130>/Constant' */
    real_T NeckPosition[6]; /* '<S118>/Neck Position' */
    real_T Constant_a[3]; /* '<S127>/Constant' */
    real_T SFunction_o[16]; /* '<S132>/S-Function' */
    real_T SFunction_a[16]; /* '<S121>/S-Function' */
    real_T Constant_j; /* '<S138>/Constant' */
    real_T SFunction_c2[16]; /* '<S133>/S-Function' */
    real_T Constant_g; /* '<S140>/Constant' */
    real_T NeckPosition_m[6]; /* '<S120>/Neck Position' */
    real_T Constant_aq[3]; /* '<S137>/Constant' */
    real_T Switch[16]; /* '<S109>/Switch' */
    real_T Sum[23]; /* '<S109>/Sum' */
    real_T SFunction_e[174]; /* '<S114>/S-Function' */
    real_T SFunction_cf[174]; /* '<S115>/S-Function' */
    real_T SFunction_on[6]; /* '<S112>/S-Function' */
    real_T SFunction_oe[16]; /* '<S25>/S-Function' */
    real_T SFunction_j[16]; /* '<S26>/S-Function' */
    real_T SFunction_c4[174]; /* '<S23>/S-Function' */
    real_T SFunction_b3[174]; /* '<S24>/S-Function' */
    real_T SFunction_ej[6]; /* '<S20>/S-Function' */
    real_T SFunction_h[6]; /* '<S21>/S-Function' */
    real_T SFunction_dk[16]; /* '<S28>/S-Function' */
    real_T SFunction_n[174]; /* '<S22>/S-Function' */
    real_T MinimumJerkTrajectoryGenerator2[3]; /* '<S33>/Minimum Jerk Trajectory Generator2' */
    real_T MinimumJerkTrajectoryGenerato_d[3]; /* '<S33>/Minimum Jerk Trajectory Generator2' */
    real_T MinimumJerkTrajectoryGenerato_b[3]; /* '<S33>/Minimum Jerk Trajectory Generator2' */
    real_T Switch_c[23]; /* '<S111>/Switch' */
    real_T Saturation[23]; /* '<Root>/Saturation' */
    real_T SFunction_o1[23]; /* '<S152>/S-Function' */
    real_T SFunction_o2[23]; /* '<S152>/S-Function' */
    real_T index1; /* '<S7>/index1' */
    real_T inRange; /* '<S7>/MATLAB Function' */
    real_T reflectedInertia[529]; /* '<S111>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}' */
    real_T QPTwoFeet_o1[12]; /* '<S146>/QP Two Feet' */
    real_T QPTwoFeet_o2; /* '<S146>/QP Two Feet' */
    real_T QPTwoFeet_o1_j[6]; /* '<S145>/QP Two Feet' */
    real_T QPTwoFeet_o2_h; /* '<S145>/QP Two Feet' */
    real_T nu_b_equivalent[6]; /* '<S109>/References for H' */
    real_T HessianMatrixQP1Foot[36]; /* '<S6>/Balancing Controller ' */
    real_T gradientQP1Foot[6]; /* '<S6>/Balancing Controller ' */
    real_T ConstraintsMatrixQP1Foot[114]; /* '<S6>/Balancing Controller ' */
    real_T bVectorConstraintsQp1Foot[19]; /* '<S6>/Balancing Controller ' */
    real_T HessianMatrixQP2Feet[144]; /* '<S6>/Balancing Controller ' */
    real_T gradientQP2Feet[12]; /* '<S6>/Balancing Controller ' */
    real_T ConstraintsMatrixQP2Feet[456]; /* '<S6>/Balancing Controller ' */
    real_T bVectorConstraintsQp2Feet[38]; /* '<S6>/Balancing Controller ' */
    real_T wR_WBDT[6]; /* '<S34>/right_foot_wrench' */
    real_T wL_WBDT[6]; /* '<S34>/left_foot_wrench' */
    real_T Constant7_m[16]; /* '<S66>/Constant7' */
    real_T jointAngles[23]; /* '<S34>/jointAngles' */
    real_T SFunction_k[16]; /* '<S82>/S-Function' */
    real_T SFunction_ao[16]; /* '<S79>/S-Function' */
    real_T Constant_f; /* '<S88>/Constant' */
    real_T SFunction_nj[16]; /* '<S83>/S-Function' */
    real_T inertial_n[12]; /* '<S34>/inertial' */
    real_T Constant_df; /* '<S90>/Constant' */
    real_T NeckPosition_p[6]; /* '<S78>/Neck Position' */
    real_T Constant_ao[3]; /* '<S87>/Constant' */
    real_T Switch6[16]; /* '<S78>/Switch6' */
    real_T SFunction_jj[16]; /* '<S106>/S-Function' */
    real_T Constant_k; /* '<S102>/Constant' */
    real_T Constant_n; /* '<S104>/Constant' */
    real_T SFunction_oz[16]; /* '<S92>/S-Function' */
    real_T SFunction_i[16]; /* '<S81>/S-Function' */
    real_T Constant_i; /* '<S98>/Constant' */
    real_T SFunction_ef[16]; /* '<S93>/S-Function' */
    real_T Constant_n1; /* '<S100>/Constant' */
    real_T NeckPosition_k[6]; /* '<S80>/Neck Position' */
    real_T Constant_p[3]; /* '<S97>/Constant' */
    real_T Switch6_g[16]; /* '<S80>/Switch6' */
    real_T SFunction_dn[16]; /* '<S107>/S-Function' */
    real_T SFunction_e0[174]; /* '<S76>/S-Function' */
    real_T SFunction_da[174]; /* '<S77>/S-Function' */
    real_T SFunction_nr[23]; /* '<S75>/S-Function' */
    real_T Constant1_e[6]; /* '<S34>/Constant1' */
    real_T TmpSignalConversionAtMinimumJer[29];
    real_T MinimumJerkTrajectoryGenerator[29]; /* '<S34>/Minimum Jerk Trajectory Generator' */
    real_T Reshape1[16]; /* '<S34>/Reshape1' */
    real_T w_H_b[16]; /* '<S34>/stateMachineYogaFCN' */
    real_T CoM_des[3]; /* '<S34>/stateMachineYogaFCN' */
    real_T qj_des[23]; /* '<S34>/stateMachineYogaFCN' */
    real_T constraints[2]; /* '<S34>/stateMachineYogaFCN' */
    real_T currentState; /* '<S34>/stateMachineYogaFCN' */
    real_T jointsSmoothingTime; /* '<S34>/stateMachineYogaFCN' */
    real_T nu_b[6]; /* '<S65>/Compute Base Velocity' */
    real_T torso[23]; /* '<S64>/DoFs converter' */
    real_T left_arm[23]; /* '<S64>/DoFs converter' */
    real_T right_arm[23]; /* '<S64>/DoFs converter' */
    real_T left_leg[23]; /* '<S64>/DoFs converter' */
    real_T right_leg[23]; /* '<S64>/DoFs converter' */
    real_T Gain_i[23]; /* '<S64>/Gain' */
    real_T torso_o[23]; /* '<S64>/DoFs converter1' */
    real_T left_arm_p[23]; /* '<S64>/DoFs converter1' */
    real_T right_arm_o[23]; /* '<S64>/DoFs converter1' */
    real_T left_leg_m[23]; /* '<S64>/DoFs converter1' */
    real_T right_leg_o[23]; /* '<S64>/DoFs converter1' */
    real_T Constant7_p[16]; /* '<S35>/Constant7' */
    real_T jointAngles_p[23]; /* '<S32>/jointAngles' */
    real_T SFunction_jo[16]; /* '<S44>/S-Function' */
    real_T SFunction_l[16]; /* '<S45>/S-Function' */
    real_T SFunction_k0[16]; /* '<S42>/S-Function' */
    real_T SFunction_bdr[16]; /* '<S43>/S-Function' */
    real_T NeckPosition_c[6]; /* '<S41>/Neck Position' */
    real_T Constant_e; /* '<S50>/Constant' */
    real_T IMUmeasurements[12]; /* '<S32>/IMU measurements' */
    real_T Constant_mz; /* '<S52>/Constant' */
    real_T Constant_ge[3]; /* '<S49>/Constant' */
    real_T Switch6_b[16]; /* '<S41>/Switch6' */
    real_T SFunction_p[174]; /* '<S57>/S-Function' */
    real_T SFunction_i0[174]; /* '<S58>/S-Function' */
    real_T Constant1_c[2]; /* '<S32>/Constant1' */
    real_T SFunction_m[23]; /* '<S56>/S-Function' */
    real_T Constant2[6]; /* '<S32>/Constant2' */
    real_T SFunction_ak[16]; /* '<S63>/S-Function' */
    real_T Constant_f5; /* '<S61>/Constant' */
    real_T Reshape1_p[16]; /* '<S32>/Reshape1' */
    real_T Switch1[9]; /* '<S32>/Switch1' */
    real_T Constant_b; /* '<S59>/Constant' */
    real_T Constant3[23]; /* '<S32>/Constant3' */
    real_T Constant4; /* '<S32>/Constant4' */
    real_T jointssmoothingTime; /* '<S32>/joints.smoothingTime' */
    real_T Constant5[3]; /* '<S32>/Constant5' */
    real_T Constant6[3]; /* '<S32>/Constant6' */
    real_T nu_b_m[6]; /* '<S36>/Compute Base Velocity' */
    real_T torso_n[23]; /* '<S10>/DoFs converter' */
    real_T left_arm_h[23]; /* '<S10>/DoFs converter' */
    real_T right_arm_p[23]; /* '<S10>/DoFs converter' */
    real_T left_leg_i[23]; /* '<S10>/DoFs converter' */
    real_T right_leg_m[23]; /* '<S10>/DoFs converter' */
    real_T SFunction_bm[23]; /* '<S12>/S-Function' */
    real_T torqueError[23]; /* '<S10>/Sum' */
    real_T torso_nn[23]; /* '<S10>/DoFs converter1' */
    real_T left_arm_i[23]; /* '<S10>/DoFs converter1' */
    real_T right_arm_i[23]; /* '<S10>/DoFs converter1' */
    real_T left_leg_h[23]; /* '<S10>/DoFs converter1' */
    real_T right_leg_g[23]; /* '<S10>/DoFs converter1' */
    real_T torso_b[23]; /* '<S10>/DoFs converter2' */
    real_T left_arm_b[23]; /* '<S10>/DoFs converter2' */
    real_T right_arm_f[23]; /* '<S10>/DoFs converter2' */
    real_T left_leg_f[23]; /* '<S10>/DoFs converter2' */
    real_T right_leg_j[23]; /* '<S10>/DoFs converter2' */
    real_T torso_g[23]; /* '<S10>/DoFs converter6' */
    real_T left_arm_pa[23]; /* '<S10>/DoFs converter6' */
    real_T right_arm_h[23]; /* '<S10>/DoFs converter6' */
    real_T left_leg_iu[23]; /* '<S10>/DoFs converter6' */
    real_T right_leg_h[23]; /* '<S10>/DoFs converter6' */
    real_T torqueError_a[23]; /* '<S10>/Sum4' */
    real_T torso_i[23]; /* '<S10>/DoFs converter7' */
    real_T left_arm_hk[23]; /* '<S10>/DoFs converter7' */
    real_T right_arm_e[23]; /* '<S10>/DoFs converter7' */
    real_T left_leg_b[23]; /* '<S10>/DoFs converter7' */
    real_T right_leg_p[23]; /* '<S10>/DoFs converter7' */
    real_T SFunction_ol[23]; /* '<S11>/S-Function' */
    real_T Gain_g[23]; /* '<S10>/Gain' */
    real_T torso_p[23]; /* '<S10>/DoFs converter3' */
    real_T left_arm_bo[23]; /* '<S10>/DoFs converter3' */
    real_T right_arm_l[23]; /* '<S10>/DoFs converter3' */
    real_T left_leg_p[23]; /* '<S10>/DoFs converter3' */
    real_T right_leg_a[23]; /* '<S10>/DoFs converter3' */
    real_T Gain3[23]; /* '<S10>/Gain3' */
    real_T torso_m[23]; /* '<S10>/DoFs converter4' */
    real_T left_arm_j[23]; /* '<S10>/DoFs converter4' */
    real_T right_arm_n[23]; /* '<S10>/DoFs converter4' */
    real_T left_leg_bl[23]; /* '<S10>/DoFs converter4' */
    real_T right_leg_ag[23]; /* '<S10>/DoFs converter4' */
    real_T Gain4[23]; /* '<S10>/Gain4' */
    real_T torso_og[23]; /* '<S10>/DoFs converter5' */
    real_T left_arm_e[23]; /* '<S10>/DoFs converter5' */
    real_T right_arm_il[23]; /* '<S10>/DoFs converter5' */
    real_T left_leg_c[23]; /* '<S10>/DoFs converter5' */
    real_T right_leg_d[23]; /* '<S10>/DoFs converter5' */
    uint8_T Constant2_j; /* '<S4>/Constant2' */
    boolean_T Compare; /* '<S30>/Compare' */
    boolean_T Compare_c; /* '<S29>/Compare' */
    boolean_T SMOOTH_DES_COM2; /* '<S33>/SMOOTH_DES_COM2' */
    boolean_T USE_IMU4EST_BASE1; /* '<S127>/USE_IMU4EST_BASE1' */
    boolean_T USE_IMU4EST_BASE; /* '<S118>/USE_IMU4EST_BASE' */
    boolean_T USE_IMU4EST_BASE1_f; /* '<S137>/USE_IMU4EST_BASE1' */
    boolean_T USE_IMU4EST_BASE_g; /* '<S120>/USE_IMU4EST_BASE' */
    boolean_T SMOOTH_DES_COM; /* '<S33>/SMOOTH_DES_COM' */
    boolean_T HiddenBuf_InsertedFor_OneFoot_a; /* '<S142>/ContactsTransition' */
    boolean_T not_a; /* '<S142>/not' */
    boolean_T HiddenBuf_InsertedFor_TwoFeet_a; /* '<S142>/not' */
    boolean_T u_g; /* '<S111>/ ' */
    boolean_T LogicalOperator1; /* '<S33>/Logical Operator1' */
    boolean_T LogicalOperator1_n; /* '<S2>/Logical Operator1' */
    boolean_T onOneFoot; /* '<S142>/ContactsTransition' */
    boolean_T USE_IMU4EST_BASE1_e; /* '<S87>/USE_IMU4EST_BASE1' */
    boolean_T USE_IMU4EST_BASE_d; /* '<S78>/USE_IMU4EST_BASE' */
    boolean_T USE_IMU4EST_BASE1_k; /* '<S97>/USE_IMU4EST_BASE1' */
    boolean_T USE_IMU4EST_BASE_a; /* '<S80>/USE_IMU4EST_BASE' */
    boolean_T Constant1_m; /* '<S35>/Constant1' */
    boolean_T USE_IMU4EST_BASE1_b; /* '<S49>/USE_IMU4EST_BASE1' */
    boolean_T USE_IMU4EST_BASE_c; /* '<S41>/USE_IMU4EST_BASE' */
    boolean_T Constant_ib; /* '<S32>/Constant' */
    B_MATLABFunction_torqueBala_f_T sf_MATLABFunction_l; /* '<S136>/MATLAB Function' */
    B_MATLABFunction_torqueBalanc_T sf_MATLABFunction_o; /* '<S135>/MATLAB Function' */
    B_fromImuToHomogeousTransform_T
        sf_fromImuToHomogeousTransfor_e; /* '<S120>/fromImuToHomogeousTransformFCN' */
    B_MATLABFunction_torqueBala_f_T sf_MATLABFunction_h; /* '<S126>/MATLAB Function' */
    B_MATLABFunction_torqueBalanc_T sf_MATLABFunction; /* '<S125>/MATLAB Function' */
    B_fromImuToHomogeousTransform_T
        sf_fromImuToHomogeousTransformF; /* '<S118>/fromImuToHomogeousTransformFCN' */
    B_MATLABFunction_torqueBala_j_T sf_MATLABFunction_d; /* '<S69>/MATLAB Function' */
    B_MATLABFunction_torqueBala_p_T sf_MATLABFunction_f; /* '<S68>/MATLAB Function' */
    B_MATLABFunction_torqueBala_f_T sf_MATLABFunction_jn; /* '<S96>/MATLAB Function' */
    B_MATLABFunction_torqueBalanc_T sf_MATLABFunction_ad; /* '<S95>/MATLAB Function' */
    B_fromImuToHomogeousTransform_T
        sf_fromImuToHomogeousTransfor_g; /* '<S80>/fromImuToHomogeousTransformFCN' */
    B_MATLABFunction_torqueBala_f_T sf_MATLABFunction_oa; /* '<S86>/MATLAB Function' */
    B_MATLABFunction_torqueBalanc_T sf_MATLABFunction_i; /* '<S85>/MATLAB Function' */
    B_fromImuToHomogeousTransform_T
        sf_fromImuToHomogeousTransfo_fb; /* '<S78>/fromImuToHomogeousTransformFCN' */
    B_MATLABFunction_torqueBala_p_T sf_MATLABFunction_a; /* '<S39>/MATLAB Function' */
    B_MATLABFunction_torqueBala_j_T sf_MATLABFunction_k; /* '<S38>/MATLAB Function' */
    B_MATLABFunction_torqueBala_f_T sf_MATLABFunction_l5; /* '<S48>/MATLAB Function' */
    B_MATLABFunction_torqueBalanc_T sf_MATLABFunction_os; /* '<S47>/MATLAB Function' */
    B_fromImuToHomogeousTransform_T
        sf_fromImuToHomogeousTransfor_f; /* '<S41>/fromImuToHomogeousTransformFCN' */
} B_torqueBalancingYoga_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct
{
    real_T state; /* '<S34>/stateMachineYogaFCN' */
    real_T tSwitch; /* '<S34>/stateMachineYogaFCN' */
    real_T w_H_fixedLink[16]; /* '<S34>/stateMachineYogaFCN' */
    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK; /* '<S151>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_c; /* '<S31>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator1; /* '<S33>/Minimum Jerk Trajectory Generator1' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_n; /* '<S18>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_h; /* '<S17>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_e; /* '<S16>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_b; /* '<S122>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d; /* '<S119>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_nq; /* '<S123>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } inertial_PWORK; /* '<S109>/inertial' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK; /* '<S118>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_k; /* '<S132>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_er; /* '<S121>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_hg; /* '<S133>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK_j; /* '<S120>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_dc; /* '<S114>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_a; /* '<S115>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_c4; /* '<S112>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_j; /* '<S25>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_p; /* '<S26>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_ne; /* '<S23>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_bp; /* '<S24>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d4; /* '<S20>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d1; /* '<S21>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_hl; /* '<S28>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_ay; /* '<S22>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator2; /* '<S33>/Minimum Jerk Trajectory Generator2' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_l; /* '<S5>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } RealTimeSynchronizer_PWORK; /* '<S155>/Real Time Synchronizer' */

    struct
    {
        void* blockPWork[2];
    } SimulatorSynchronizer_PWORK; /* '<S154>/Simulator Synchronizer' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_g; /* '<S152>/S-Function' */

    void* Assertion_slioAccessor; /* '<S7>/Assertion' */
    struct
    {
        void* blockPWork[2];
    } QPTwoFeet_PWORK; /* '<S146>/QP Two Feet' */

    struct
    {
        void* blockPWork[2];
    } QPTwoFeet_PWORK_l; /* '<S145>/QP Two Feet' */

    struct
    {
        void* blockPWork[2];
    } right_foot_wrench_PWORK; /* '<S34>/right_foot_wrench' */

    struct
    {
        void* blockPWork[2];
    } left_foot_wrench_PWORK; /* '<S34>/left_foot_wrench' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_pz; /* '<S82>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_bg; /* '<S79>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_m; /* '<S83>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } inertial_PWORK_a; /* '<S34>/inertial' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK_n; /* '<S78>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_hk; /* '<S106>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_kk; /* '<S92>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_lc; /* '<S81>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_m2; /* '<S93>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK_p; /* '<S80>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_d0; /* '<S107>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_j2; /* '<S76>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_jg; /* '<S77>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_bs; /* '<S75>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } MinimumJerkTrajectoryGenerator_; /* '<S34>/Minimum Jerk Trajectory Generator' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter_PWORK; /* '<S64>/DoFs converter' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter1_PWORK; /* '<S64>/DoFs converter1' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_ga; /* '<S44>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_es; /* '<S45>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_o; /* '<S42>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_kj; /* '<S43>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } NeckPosition_PWORK_c; /* '<S41>/Neck Position' */

    struct
    {
        void* blockPWork[2];
    } IMUmeasurements_PWORK; /* '<S32>/IMU measurements' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_l5; /* '<S57>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_i; /* '<S58>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_j1; /* '<S56>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_gx; /* '<S63>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter_PWORK_o; /* '<S10>/DoFs converter' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_f; /* '<S12>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter1_PWORK_j; /* '<S10>/DoFs converter1' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter2_PWORK; /* '<S10>/DoFs converter2' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter6_PWORK; /* '<S10>/DoFs converter6' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter7_PWORK; /* '<S10>/DoFs converter7' */

    struct
    {
        void* blockPWork[2];
    } SFunction_PWORK_dk; /* '<S11>/S-Function' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter3_PWORK; /* '<S10>/DoFs converter3' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter4_PWORK; /* '<S10>/DoFs converter4' */

    struct
    {
        void* blockPWork[2];
    } DoFsconverter5_PWORK; /* '<S10>/DoFs converter5' */

    boolean_T state_not_empty; /* '<S34>/stateMachineYogaFCN' */
    boolean_T tSwitch_not_empty; /* '<S34>/stateMachineYogaFCN' */
    boolean_T w_H_fixedLink_not_empty; /* '<S34>/stateMachineYogaFCN' */
    boolean_T emergencystopjointlimits_MODE; /* '<Root>/emergency stop: joint limits' */
    boolean_T TwoFeet_MODE; /* '<S142>/Two Feet' */
    boolean_T OneFoot_MODE; /* '<S142>/One Foot' */
    boolean_T StateMachineYoga_MODE; /* '<S4>/State Machine Yoga' */
    boolean_T VisualizeGainTuning_MODE; /* '<S33>/Visualize Gain Tuning  ' */
    boolean_T InternalCoordinator_MODE; /* '<S4>/Internal Coordinator' */
    boolean_T Visualizer_MODE; /* '<S2>/Visualizer' */
    DW_MATLABFunction_torqueBal_k_T sf_MATLABFunction_l; /* '<S136>/MATLAB Function' */
    DW_MATLABFunction_torqueBalan_T sf_MATLABFunction_o; /* '<S135>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_k_T sf_MATLABFunction_h; /* '<S126>/MATLAB Function' */
    DW_MATLABFunction_torqueBalan_T sf_MATLABFunction; /* '<S125>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_l_T sf_MATLABFunction_d; /* '<S69>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_e_T sf_MATLABFunction_f; /* '<S68>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_k_T sf_MATLABFunction_jn; /* '<S96>/MATLAB Function' */
    DW_MATLABFunction_torqueBalan_T sf_MATLABFunction_ad; /* '<S95>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_k_T sf_MATLABFunction_oa; /* '<S86>/MATLAB Function' */
    DW_MATLABFunction_torqueBalan_T sf_MATLABFunction_i; /* '<S85>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_e_T sf_MATLABFunction_a; /* '<S39>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_l_T sf_MATLABFunction_k; /* '<S38>/MATLAB Function' */
    DW_MATLABFunction_torqueBal_k_T sf_MATLABFunction_l5; /* '<S48>/MATLAB Function' */
    DW_MATLABFunction_torqueBalan_T sf_MATLABFunction_os; /* '<S47>/MATLAB Function' */
} DW_torqueBalancingYoga_T;

/* Parameters (auto storage) */
struct P_torqueBalancingYoga_T_
{
    struct_xTYJupe8umqLUY4vWPaUnH Sm; /* Variable: Sm
                                       * Referenced by:
                                       *   '<S32>/joints.smoothingTime'
                                       *   '<S34>/stateMachineYogaFCN'
                                       */
    struct_nKwr8yfDbWZ5195E5duazE Config; /* Variable: Config
                                           * Referenced by:
                                           *   '<Root>/ON_GAZEBO 1'
                                           *   '<S2>/ON_GAZEBO 3'
                                           *   '<S2>/ON_GAZEBO 4'
                                           *   '<S8>/ON_GAZEBO '
                                           *   '<S32>/Reference Generator CoM'
                                           *   '<S32>/Constant'
                                           *   '<S32>/Constant1'
                                           *   '<S33>/ON_GAZEBO 3'
                                           *   '<S33>/ON_GAZEBO 4'
                                           *   '<S33>/SMOOTH_DES_COM'
                                           *   '<S33>/SMOOTH_DES_COM2'
                                           *   '<S34>/ON_GAZEBO 3'
                                           *   '<S34>/ON_GAZEBO 4'
                                           *   '<S111>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}'
                                           *   '<S111>/ '
                                           *   '<S111>/Gain'
                                           *   '<S15>/Add motor reflected inertias'
                                           *   '<S35>/Constant1'
                                           *   '<S41>/fromImuToHomogeousTransformFCN'
                                           *   '<S41>/USE_IMU4EST_BASE'
                                           *   '<S78>/fromImuToHomogeousTransformFCN'
                                           *   '<S78>/USE_IMU4EST_BASE'
                                           *   '<S80>/fromImuToHomogeousTransformFCN'
                                           *   '<S80>/USE_IMU4EST_BASE'
                                           *   '<S118>/fromImuToHomogeousTransformFCN'
                                           *   '<S118>/USE_IMU4EST_BASE'
                                           *   '<S120>/fromImuToHomogeousTransformFCN'
                                           *   '<S120>/USE_IMU4EST_BASE'
                                           *   '<S49>/USE_IMU4EST_BASE1'
                                           *   '<S87>/USE_IMU4EST_BASE1'
                                           *   '<S97>/USE_IMU4EST_BASE1'
                                           *   '<S127>/USE_IMU4EST_BASE1'
                                           *   '<S137>/USE_IMU4EST_BASE1'
                                           */
    struct_r750FiCaocHqTuxb1nF3DC Gain; /* Variable: Gain
                                         * Referenced by:
                                         *   '<S6>/Balancing Controller '
                                         *   '<S34>/stateMachineYogaFCN'
                                         */
    struct_iWQaDqPVqHoUwXRopcYWkG Reg; /* Variable: Reg
                                        * Referenced by:
                                        *   '<S6>/Balancing Controller '
                                        *   '<S109>/References for H'
                                        *   '<S36>/Compute Base Velocity'
                                        *   '<S65>/Compute Base Velocity'
                                        */
    real_T ConstraintsMatrix[114]; /* Variable: ConstraintsMatrix
                                    * Referenced by: '<S6>/Constant'
                                    */
    real_T ROBOT_DOF_FOR_SIMULINK[529]; /* Variable: ROBOT_DOF_FOR_SIMULINK
                                         * Referenced by: '<S6>/    5'
                                         */
    real_T bVectorConstraints[19]; /* Variable: bVectorConstraints
                                    * Referenced by: '<S6>/Constant1'
                                    */
    struct_5rbQfDr6lxmxJDTJAmE7vH Sat; /* Variable: Sat
                                        * Referenced by: '<Root>/Saturation'
                                        */
    real_T CompareToConstant_const; /* Mask Parameter: CompareToConstant_const
                                     * Referenced by: '<S50>/Constant'
                                     */
    real_T CompareToConstant_const_j; /* Mask Parameter: CompareToConstant_const_j
                                       * Referenced by: '<S52>/Constant'
                                       */
    real_T CompareToConstant_const_o; /* Mask Parameter: CompareToConstant_const_o
                                       * Referenced by: '<S61>/Constant'
                                       */
    real_T CompareToConstant_const_ok; /* Mask Parameter: CompareToConstant_const_ok
                                        * Referenced by: '<S59>/Constant'
                                        */
    real_T CompareToConstant_const_p; /* Mask Parameter: CompareToConstant_const_p
                                       * Referenced by: '<S88>/Constant'
                                       */
    real_T CompareToConstant_const_od; /* Mask Parameter: CompareToConstant_const_od
                                        * Referenced by: '<S90>/Constant'
                                        */
    real_T CompareToConstant_const_ph; /* Mask Parameter: CompareToConstant_const_ph
                                        * Referenced by: '<S102>/Constant'
                                        */
    real_T CompareToConstant_const_pf; /* Mask Parameter: CompareToConstant_const_pf
                                        * Referenced by: '<S104>/Constant'
                                        */
    real_T CompareToConstant_const_f; /* Mask Parameter: CompareToConstant_const_f
                                       * Referenced by: '<S98>/Constant'
                                       */
    real_T CompareToConstant_const_e; /* Mask Parameter: CompareToConstant_const_e
                                       * Referenced by: '<S100>/Constant'
                                       */
    real_T CompareToConstant_const_h; /* Mask Parameter: CompareToConstant_const_h
                                       * Referenced by: '<S128>/Constant'
                                       */
    real_T CompareToConstant_const_c; /* Mask Parameter: CompareToConstant_const_c
                                       * Referenced by: '<S130>/Constant'
                                       */
    real_T CompareToConstant_const_ea; /* Mask Parameter: CompareToConstant_const_ea
                                        * Referenced by: '<S138>/Constant'
                                        */
    real_T CompareToConstant_const_e5; /* Mask Parameter: CompareToConstant_const_e5
                                        * Referenced by: '<S140>/Constant'
                                        */
    uint8_T Coordinator_BitMask; /* Mask Parameter: Coordinator_BitMask
                                  * Referenced by: '<S4>/Coordinator'
                                  */
    uint8_T Yoga_BitMask; /* Mask Parameter: Yoga_BitMask
                           * Referenced by: '<S4>/Yoga'
                           */
    real_T Gain_Gain; /* Expression: 180/pi
                       * Referenced by: '<S10>/Gain'
                       */
    real_T Gain3_Gain; /* Expression: 180/pi
                        * Referenced by: '<S10>/Gain3'
                        */
    real_T Gain4_Gain; /* Expression: 180/pi
                        * Referenced by: '<S10>/Gain4'
                        */
    real_T Gain_Gain_h; /* Expression: pi/180
                         * Referenced by: '<S49>/Gain'
                         */
    real_T Constant7_Value[16]; /* Expression: eye(4)
                                 * Referenced by: '<S35>/Constant7'
                                 */
    real_T Constant_Value[3]; /* Expression: zeros(3,1)
                               * Referenced by: '<S49>/Constant'
                               */
    real_T Constant2_Value[6]; /* Expression: zeros(6,1)
                                * Referenced by: '<S32>/Constant2'
                                */
    real_T Constant3_Value[23]; /* Expression: Gain.impedances(1,:)
                                 * Referenced by: '<S32>/Constant3'
                                 */
    real_T Constant4_Value; /* Expression: 1
                             * Referenced by: '<S32>/Constant4'
                             */
    real_T Constant5_Value[3]; /* Expression: diag(Gain.KP_COM)
                                * Referenced by: '<S32>/Constant5'
                                */
    real_T Constant6_Value[3]; /* Expression: diag(Gain.KD_COM)
                                * Referenced by: '<S32>/Constant6'
                                */
    real_T Gain_Gain_e; /* Expression: 180/pi
                         * Referenced by: '<S64>/Gain'
                         */
    real_T Gain_Gain_a; /* Expression: pi/180
                         * Referenced by: '<S87>/Gain'
                         */
    real_T Gain_Gain_hh; /* Expression: pi/180
                          * Referenced by: '<S97>/Gain'
                          */
    real_T Constant7_Value_b[16]; /* Expression: eye(4)
                                   * Referenced by: '<S66>/Constant7'
                                   */
    real_T Constant_Value_c[3]; /* Expression: zeros(3,1)
                                 * Referenced by: '<S87>/Constant'
                                 */
    real_T Constant_Value_cj[3]; /* Expression: zeros(3,1)
                                  * Referenced by: '<S97>/Constant'
                                  */
    real_T Constant1_Value[6]; /* Expression: zeros(6,1)
                                * Referenced by: '<S34>/Constant1'
                                */
    real_T Gain_Gain_hk; /* Expression: pi/180
                          * Referenced by: '<S127>/Gain'
                          */
    real_T Gain_Gain_c; /* Expression: pi/180
                         * Referenced by: '<S137>/Gain'
                         */
    real_T index1_Value; /* Expression: 0.01
                          * Referenced by: '<S7>/index1'
                          */
    real_T Constant_Value_d[6]; /* Expression: zeros(6,1)
                                 * Referenced by: '<S17>/Constant'
                                 */
    real_T Gain_Gain_o; /* Expression: 0
                         * Referenced by: '<S17>/Gain'
                         */
    real_T Constant7_Value_o[16]; /* Expression: eye(4)
                                   * Referenced by: '<S113>/Constant7'
                                   */
    real_T Constant_Value_o[3]; /* Expression: zeros(3,1)
                                 * Referenced by: '<S127>/Constant'
                                 */
    real_T Constant_Value_i[3]; /* Expression: zeros(3,1)
                                 * Referenced by: '<S137>/Constant'
                                 */
    real_T Switch_Threshold; /* Expression: 0.1
                              * Referenced by: '<S109>/Switch'
                              */
    uint8_T Constant_Value_n; /* Computed Parameter: Constant_Value_n
                               * Referenced by: '<S30>/Constant'
                               */
    uint8_T Constant_Value_it; /* Expression: Sm.SM_TYPE_BIN
                                * Referenced by: '<S4>/Constant'
                                */
    uint8_T Constant_Value_dx; /* Computed Parameter: Constant_Value_dx
                                * Referenced by: '<S29>/Constant'
                                */
    uint8_T Constant2_Value_o; /* Expression: Sm.SM_TYPE_BIN
                                * Referenced by: '<S4>/Constant2'
                                */
};

/* Real-time Model Data Structure */
struct tag_RTM_torqueBalancingYoga_T
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

/* Class declaration for model torqueBalancingYoga */
class torqueBalancingYogaModelClass
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
    torqueBalancingYogaModelClass();

    /* Destructor */
    ~torqueBalancingYogaModelClass();

    /* Real-Time Model get method */
    RT_MODEL_torqueBalancingYoga_T* getRTM();

    /* private data and function members */
private:
    /* Tunable parameters */
    P_torqueBalancingYoga_T torqueBalancingYoga_P;

    /* Block signals */
    B_torqueBalancingYoga_T torqueBalancingYoga_B;

    /* Block states */
    DW_torqueBalancingYoga_T torqueBalancingYoga_DW;

    /* Real-Time Model */
    RT_MODEL_torqueBalancingYoga_T torqueBalancingYoga_M;

    /* private member function(s) for subsystem '<Root>'*/
    void torqueBalancingYoga_xswap_g(real_T x[529], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xgetrf_m(real_T A[529], int32_T ipiv[23], int32_T* info);
    void torqueBalancingYoga_xtrsm_l(const real_T A[529], real_T B[529]);
    void torqueBalancingYoga_xtrsm_lc(const real_T A[529], real_T B[529]);
    void torqueBalancingYoga_mrdivide_g(const real_T A[529], const real_T B[529], real_T y[529]);
    void torqueBalancingYoga_eye(real_T I[16]);
    void torqueBalancingYoga_xgetrf_f(real_T A[16], int32_T ipiv[4], int32_T* info);
    void torqueBalancingYoga_mrdivide_l(real_T A[16], const real_T B[16]);
    void torqueBalancingYoga_mldivide_a(const real_T A[16], real_T B[4]);
    real_T torqueBalancingYoga_norm(const real_T x[6]);
    real_T torqueBalancingYoga_norm_p(const real_T x[2]);
    void torqueBalancingYoga_mldivide_m(const real_T A[36], real_T B[72]);
    void torqueBalancingYoga_mldivide(const real_T A[36], real_T B[72]);
    void torqueBalancingYoga_xswap(real_T x[529], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xgetrf_l(real_T A[529], int32_T ipiv[23], int32_T* info);
    void torqueBalancingYoga_xtrsm(const real_T A[529], real_T B[529]);
    void torqueBalancingYoga_xtrsm_d(const real_T A[529], real_T B[529]);
    void torqueBalancingYoga_mrdivide(const real_T A[529], const real_T B[529], real_T y[529]);
    void torqueBala_computeMotorsInertia(const struct_nKwr8yfDbWZ5195E5duazE* b_Config,
                                         real_T reflectedInertia[529]);
    void torqueBalancingYoga_mldivide_e(const real_T A[36], real_T B[36]);
    void torqueBalancingYoga_xzgetrf(real_T A[36], int32_T ipiv[6], int32_T* info);
    void torqueBalancingYoga_inv(const real_T x[36], real_T y[36]);
    real_T torqueBalancingYoga_xnrm2(int32_T n, const real_T x[72], int32_T ix0);
    real_T torqueBalancingYoga_xnrm2_f(int32_T n, const real_T x[6], int32_T ix0);
    void torqueBalancingYoga_xaxpy_kt(int32_T n,
                                      real_T a,
                                      const real_T x[12],
                                      int32_T ix0,
                                      real_T y[72],
                                      int32_T iy0);
    void torqueBalancingYoga_xaxpy_k(int32_T n,
                                     real_T a,
                                     const real_T x[72],
                                     int32_T ix0,
                                     real_T y[12],
                                     int32_T iy0);
    real_T torqueBalancingYoga_xdotc(int32_T n,
                                     const real_T x[72],
                                     int32_T ix0,
                                     const real_T y[72],
                                     int32_T iy0);
    void torqueBalancingYoga_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[72], int32_T iy0);
    real_T torqueBalancingYoga_xdotc_p(int32_T n,
                                       const real_T x[36],
                                       int32_T ix0,
                                       const real_T y[36],
                                       int32_T iy0);
    void torqueBalancingYoga_xaxpy_kto(int32_T n, real_T a, int32_T ix0, real_T y[36], int32_T iy0);
    void torqueBalancingYoga_xscal(real_T a, real_T x[72], int32_T ix0);
    void torqueBalancingYoga_xscal_k(real_T a, real_T x[36], int32_T ix0);
    void torqueBalancingYoga_xswap_m(real_T x[36], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xswap_mk(real_T x[72], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xrotg(real_T* a, real_T* b, real_T* c, real_T* s);
    void torqueBalancingYoga_xrot(real_T x[36], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void torqueBalancingYoga_xrot_c(real_T x[72], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void torqueBalancingYoga_svd(const real_T A[72], real_T U[72], real_T s[6], real_T V[36]);
    void torqueBalancingYoga_pinv(const real_T A[72], real_T tol, real_T X[72]);
    void torqueBalancingYoga_eye_m(real_T I[144]);
    void torqueBalancingYoga_xswap_mks(real_T x[841], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xgetrf_g(real_T A[841], int32_T ipiv[29], int32_T* info);
    void torqueBalancingYoga_xtrsm_b(const real_T A[841], real_T B[348]);
    void torqueBalancingYoga_xtrsm_bp(const real_T A[841], real_T B[348]);
    void torqueBalancingYoga_mrdivide_a(const real_T A[348], const real_T B[841], real_T y[348]);
    void torqueBalancingYoga_mrdivide_al(real_T A[138], const real_T B[36]);
    void torqueBalancingYoga_pinvDamped(const real_T A[276], real_T regDamp, real_T pinvDampA[276]);
    void torqueBalancingYoga_eye_mj(real_T I[529]);
    void torqueBalancingYoga_blkdiag(const real_T varargin_1[9],
                                     const real_T varargin_2[9],
                                     real_T y[36]);
    void torqueBalancingYoga_blkdiag_l(const real_T varargin_1[114],
                                       const real_T varargin_2[114],
                                       real_T y[456]);
    void torqueBalancingYoga_diag(const real_T v[23], real_T d[529]);
    real_T torqueBalancingYoga_xnrm2_fe(int32_T n, const real_T x[529], int32_T ix0);
    real_T torqueBalancingYoga_xnrm2_fes(int32_T n, const real_T x[23], int32_T ix0);
    void torqueBalancingYog_xaxpy_ktoigf(int32_T n,
                                         real_T a,
                                         const real_T x[23],
                                         int32_T ix0,
                                         real_T y[529],
                                         int32_T iy0);
    void torqueBalancingYoga_xaxpy_ktoig(int32_T n,
                                         real_T a,
                                         const real_T x[529],
                                         int32_T ix0,
                                         real_T y[23],
                                         int32_T iy0);
    real_T torqueBalancingYoga_xdotc_pe(int32_T n,
                                        const real_T x[529],
                                        int32_T ix0,
                                        const real_T y[529],
                                        int32_T iy0);
    void
    torqueBalancingYoga_xaxpy_ktoi(int32_T n, real_T a, int32_T ix0, real_T y[529], int32_T iy0);
    void torqueBalancingYoga_xscal_kg(real_T a, real_T x[529], int32_T ix0);
    void torqueBalancingYoga_xswap_mksh(real_T x[529], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xrot_cm(real_T x[529], int32_T ix0, int32_T iy0, real_T c, real_T s);
    void torqueBalancingYoga_svd_g(const real_T A[529], real_T U[529], real_T s[23], real_T V[529]);
    void
    torqueBalancingYoga_xgemm(int32_T k, const real_T A[529], const real_T B[529], real_T C[529]);
    void torqueBalancingYoga_pinv_d(const real_T A[529], real_T tol, real_T X[529]);
    void torqueBalancingYoga_diag_h(const real_T v[23], real_T d[529]);
    void torqueBalancingYoga_xswap_mkshm(real_T x[529], int32_T ix0, int32_T iy0);
    void torqueBalancingYoga_xgetrf_gr(real_T A[529], int32_T ipiv[23], int32_T* info);
    void torqueBalancingYoga_xtrsm_bps(const real_T A[529], real_T B[276]);
    void torqueBalancingYoga_xtrsm_bpsu(const real_T A[529], real_T B[276]);
    void
    torqueBalancingYog_pinvDamped_l(const real_T A[276], real_T regDamp, real_T pinvDampA[276]);
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
 * '<Root>' : 'torqueBalancingYoga'
 * '<S1>'   : 'torqueBalancingYoga/Configuration'
 * '<S2>'   : 'torqueBalancingYoga/Dump and visualize'
 * '<S3>'   : 'torqueBalancingYoga/Dynamics and Kinematics'
 * '<S4>'   : 'torqueBalancingYoga/Robot State and References'
 * '<S5>'   : 'torqueBalancingYoga/Set References'
 * '<S6>'   : 'torqueBalancingYoga/controller_QP'
 * '<S7>'   : 'torqueBalancingYoga/emergency stop: joint limits'
 * '<S8>'   : 'torqueBalancingYoga/synchronizer'
 * '<S9>'   : 'torqueBalancingYoga/tauDot Saturation'
 * '<S10>'  : 'torqueBalancingYoga/Dump and visualize/Visualizer'
 * '<S11>'  : 'torqueBalancingYoga/Dump and visualize/Visualizer/Get Measurement'
 * '<S12>'  : 'torqueBalancingYoga/Dump and visualize/Visualizer/Get Measurement1'
 * '<S13>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics'
 * '<S14>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics'
 * '<S15>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Add motors reflected inertias'
 * '<S16>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Centroidal Momentum'
 * '<S17>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Get Bias Forces'
 * '<S18>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Mass Matrix'
 * '<S19>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Add motors reflected
 * inertias/Add motor reflected inertias'
 * '<S20>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/DotJ Nu l_sole '
 * '<S21>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/DotJ Nu r_sole   '
 * '<S22>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian com'
 * '<S23>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian l_sole'
 * '<S24>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian r_sole'
 * '<S25>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/l_sole'
 * '<S26>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/r_sole'
 * '<S27>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/xCom'
 * '<S28>'  : 'torqueBalancingYoga/Dynamics and Kinematics/Kinematics/xCom/CoM'
 * '<S29>'  : 'torqueBalancingYoga/Robot State and References/Compare To Zero'
 * '<S30>'  : 'torqueBalancingYoga/Robot State and References/Compare To Zero1'
 * '<S31>'  : 'torqueBalancingYoga/Robot State and References/Get Measurement'
 * '<S32>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator'
 * '<S33>'  : 'torqueBalancingYoga/Robot State and References/Select State and References'
 * '<S34>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga'
 * '<S35>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link'
 * '<S36>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity'
 * '<S37>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Reference
 * Generator CoM'
 * '<S38>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder '
 * '<S39>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder 1'
 * '<S40>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/xCom'
 * '<S41>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform '
 * '<S42>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to world transform (fixed base)'
 * '<S43>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/RFoot to world transform (fixed base)'
 * '<S44>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /Fixed base to imu transform'
 * '<S45>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /Fixed base to root link transform'
 * '<S46>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /fromImuToHomogeousTransformFCN'
 * '<S47>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 1'
 * '<S48>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 2'
 * '<S49>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /neck transform'
 * '<S50>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 1/Compare To Constant'
 * '<S51>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 1/MATLAB Function'
 * '<S52>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 2/Compare To Constant'
 * '<S53>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to
 * fixed_link/LFoot to base link transform /holder 2/MATLAB Function'
 * '<S54>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity/Compute Base Velocity'
 * '<S55>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity/Feet Jacobians'
 * '<S56>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity/Get Measurement'
 * '<S57>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity/Feet Jacobians/Jacobian LFoot'
 * '<S58>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute State
 * Velocity/Feet Jacobians/Jacobian RFoot'
 * '<S59>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder /Compare
 * To Constant'
 * '<S60>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder /MATLAB
 * Function'
 * '<S61>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder 1/Compare
 * To Constant'
 * '<S62>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/holder 1/MATLAB
 * Function'
 * '<S63>'  : 'torqueBalancingYoga/Robot State and References/Internal Coordinator/xCom/CoM'
 * '<S64>'  : 'torqueBalancingYoga/Robot State and References/Select State and References/Visualize
 * Gain Tuning  '
 * '<S65>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity'
 * '<S66>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform'
 * '<S67>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Visualise extrnal
 * wrenches  '
 * '<S68>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 1'
 * '<S69>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 2'
 * '<S70>'  : 'torqueBalancingYoga/Robot State and References/State Machine
 * Yoga/stateMachineYogaFCN'
 * '<S71>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom'
 * '<S72>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom1'
 * '<S73>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity/Compute Base Velocity'
 * '<S74>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity/Feet Jacobians'
 * '<S75>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity/Get Measurement'
 * '<S76>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity/Feet Jacobians/Jacobian LFoot'
 * '<S77>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute State
 * Velocity/Feet Jacobians/Jacobian RFoot'
 * '<S78>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform '
 * '<S79>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to world transform (fixed base)'
 * '<S80>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform'
 * '<S81>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to world transform (fixed base)'
 * '<S82>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /Fixed base to imu transform'
 * '<S83>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /Fixed base to root link transform'
 * '<S84>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /fromImuToHomogeousTransformFCN'
 * '<S85>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1'
 * '<S86>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2'
 * '<S87>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /neck transform'
 * '<S88>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1/Compare To Constant'
 * '<S89>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1/MATLAB Function'
 * '<S90>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2/Compare To Constant'
 * '<S91>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2/MATLAB Function'
 * '<S92>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/Fixed base to imu transform'
 * '<S93>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/Fixed base to root link transform'
 * '<S94>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/fromImuToHomogeousTransformFCN'
 * '<S95>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1'
 * '<S96>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2'
 * '<S97>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/neck transform'
 * '<S98>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1/Compare To Constant'
 * '<S99>'  : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1/MATLAB Function'
 * '<S100>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2/Compare To Constant'
 * '<S101>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2/MATLAB Function'
 * '<S102>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 1/Compare To
 * Constant'
 * '<S103>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 1/MATLAB
 * Function'
 * '<S104>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 2/Compare To
 * Constant'
 * '<S105>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/holder 2/MATLAB
 * Function'
 * '<S106>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom/CoM'
 * '<S107>' : 'torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom1/CoM'
 * '<S108>' : 'torqueBalancingYoga/controller_QP/Balancing Controller '
 * '<S109>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral'
 * '<S110>' : 'torqueBalancingYoga/controller_QP/Compute joint torques'
 * '<S111>' : 'torqueBalancingYoga/controller_QP/Compute joint torques (motor reflected inertia)'
 * '<S112>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Centroidal
 * Momentum'
 * '<S113>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform'
 * '<S114>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Jacobian'
 * '<S115>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Jacobian1'
 * '<S116>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/References for H'
 * '<S117>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/choose base to
 * world transform'
 * '<S118>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform '
 * '<S119>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to world transform (fixed base)'
 * '<S120>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform'
 * '<S121>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to world transform (fixed base)'
 * '<S122>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /Fixed base to imu transform'
 * '<S123>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /Fixed base to root link transform'
 * '<S124>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /fromImuToHomogeousTransformFCN'
 * '<S125>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1'
 * '<S126>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2'
 * '<S127>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /neck transform'
 * '<S128>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1/Compare To Constant'
 * '<S129>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 1/MATLAB Function'
 * '<S130>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2/Compare To Constant'
 * '<S131>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/LFoot to base link transform /holder 2/MATLAB Function'
 * '<S132>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/Fixed base to imu transform'
 * '<S133>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/Fixed base to root link transform'
 * '<S134>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/fromImuToHomogeousTransformFCN'
 * '<S135>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1'
 * '<S136>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2'
 * '<S137>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/neck transform'
 * '<S138>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1/Compare To Constant'
 * '<S139>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 1/MATLAB Function'
 * '<S140>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2/Compare To Constant'
 * '<S141>' : 'torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to
 * fixed link transform/RFoot to base link transform/holder 2/MATLAB Function'
 * '<S142>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver'
 * '<S143>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/ContactsTransition'
 * '<S144>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/MATLAB Function'
 * '<S145>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot'
 * '<S146>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet'
 * '<S147>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Visualize eventual
 * QP failures'
 * '<S148>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot/Match
 * Signal Sizes1'
 * '<S149>' : 'torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet/Match
 * Signal Sizes'
 * '<S150>' : 'torqueBalancingYoga/controller_QP/Compute joint torques (motor reflected
 * inertia)/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}'
 * '<S151>' : 'torqueBalancingYoga/controller_QP/Compute joint torques (motor reflected inertia)/Get
 * Measurement'
 * '<S152>' : 'torqueBalancingYoga/emergency stop: joint limits/Get Limits'
 * '<S153>' : 'torqueBalancingYoga/emergency stop: joint limits/MATLAB Function'
 * '<S154>' : 'torqueBalancingYoga/synchronizer/GAZEBO_SYNCHRONIZER'
 * '<S155>' : 'torqueBalancingYoga/synchronizer/REAL_TIME_SYNC'
 * '<S156>' : 'torqueBalancingYoga/tauDot Saturation/Saturate the Torque Derivative'
 */
#endif /* RTW_HEADER_torqueBalancingYoga_h_ */
