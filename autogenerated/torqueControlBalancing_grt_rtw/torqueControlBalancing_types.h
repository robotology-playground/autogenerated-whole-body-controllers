/*
 * torqueControlBalancing_types.h
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

#ifndef RTW_HEADER_torqueControlBalancing_types_h_
#define RTW_HEADER_torqueControlBalancing_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#ifndef DEFINED_TYPEDEF_FOR_struct_cik6KQn5YVVa3skh1Wn00_
#define DEFINED_TYPEDEF_FOR_struct_cik6KQn5YVVa3skh1Wn00_

typedef struct
{
    real_T joints_pauseBetweenYogaMoves;
    real_T wrench_thresholdContactOn;
    real_T wrench_thresholdContactOff;
    real_T CoM_threshold;
    real_T joints_thresholdNotInContact;
    real_T joints_thresholdInContact;
    real_T initialState;
    real_T tBalancing;
    real_T tBalancingBeforeYoga;
    boolean_T yogaExtended;
    boolean_T skipYoga;
    boolean_T demoOnlyBalancing;
    boolean_T demoStartsOnRightSupport;
    boolean_T yogaAlsoOnRightFoot;
    boolean_T twoFeetYogaInLoop;
    boolean_T oneFootYogaInLoop;
    real_T yogaCounter;
    real_T CoMSmoothingTime[13];
    real_T jointsSmoothingTime[13];
    real_T scaleFactorSmoothingTime;
    real_T CoM_delta[39];
    real_T joints_references[299];
    real_T joints_leftYogaRef[624];
    real_T joints_rightYogaRef[624];
} struct_cik6KQn5YVVa3skh1Wn00;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_8fk7vJtSw4UnB4xn6cTd0B_
#define DEFINED_TYPEDEF_FOR_struct_8fk7vJtSw4UnB4xn6cTd0B_

typedef struct
{
    real_T SIMULATION_TIME;
    real_T tStep;
    boolean_T SCOPES_WRENCHES;
    boolean_T SCOPES_GAIN_SCHE_INFO;
    boolean_T SCOPES_MAIN;
    boolean_T SCOPES_QP;
    boolean_T SAVE_WORKSPACE;
    boolean_T CHECK_INTEGRATION_TIME;
    boolean_T PLOT_INTEGRATION_TIME;
    boolean_T ON_GAZEBO;
    real_T GRAV_ACC;
    real_T numOfJointsForEachControlboard[5];
    boolean_T SATURATE_TORQUE_DERIVATIVE;
    boolean_T EMERGENCY_STOP_WITH_JOINTS_LIMITS;
    boolean_T EMERGENCY_STOP_WITH_ENCODER_SPIKES;
    boolean_T USE_MOTOR_REFLECTED_INERTIA;
    boolean_T INCLUDE_COUPLING;
    boolean_T INCLUDE_HARMONIC_DRIVE_INERTIA;
    boolean_T USE_IMU4EST_BASE;
    boolean_T FILTER_IMU_YAW;
    boolean_T CORRECT_NECK_IMU;
    boolean_T USE_QP_SOLVER;
    boolean_T LEFT_RIGHT_MOVEMENTS;
    boolean_T SMOOTH_COM_DES;
    boolean_T SMOOTH_JOINT_DES;
    real_T Gamma[529];
    real_T I_m[529];
    real_T T[529];
    real_T K_ff;
    boolean_T USE_DES_JOINT_ACC_FOR_MOTORS_INERTIA;
    real_T SmoothingTimeGainScheduling;
    real_T noOscillationTime;
    real_T directionOfOscillation[3];
    real_T amplitudeOfOscillation;
    real_T frequencyOfOscillation;
    boolean_T COORDINATOR_DEMO;
} struct_8fk7vJtSw4UnB4xn6cTd0B;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_Yt0Pn6XEWtmZM9M9GtOwoD_
#define DEFINED_TYPEDEF_FOR_struct_Yt0Pn6XEWtmZM9M9GtOwoD_

typedef struct
{
    real_T pinvDamp_baseVel;
    real_T pinvDamp;
    real_T pinvTol;
    real_T KP_postural;
    real_T KD_postural;
    real_T HessianQP;
} struct_Yt0Pn6XEWtmZM9M9GtOwoD;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_glrlyr9IieN6zDechOWztB_
#define DEFINED_TYPEDEF_FOR_struct_glrlyr9IieN6zDechOWztB_

typedef struct
{
    real_T KP_CoM[39];
    real_T KD_CoM[39];
    real_T KI_AngularMomentum;
    real_T KP_AngularMomentum;
    real_T KP_postural[299];
    real_T KD_postural[23];
} struct_glrlyr9IieN6zDechOWztB;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_8ovnsQmVgsrxBHMDpqGPN_
#define DEFINED_TYPEDEF_FOR_struct_8ovnsQmVgsrxBHMDpqGPN_

typedef struct
{
    real_T torque;
    real_T uDotMax;
    real_T maxJointsPositionDelta;
} struct_8ovnsQmVgsrxBHMDpqGPN;

#endif

/* Parameters (default storage) */
typedef struct P_torqueControlBalancing_T_ P_torqueControlBalancing_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_torqueControlBalancin_T RT_MODEL_torqueControlBalanci_T;

#endif /* RTW_HEADER_torqueControlBalancing_types_h_ */
