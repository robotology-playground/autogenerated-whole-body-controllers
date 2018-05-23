/*
 * torqueBalancingYoga_types.h
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

#ifndef RTW_HEADER_torqueBalancingYoga_types_h_
#define RTW_HEADER_torqueBalancingYoga_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#ifndef DEFINED_TYPEDEF_FOR_struct_nKwr8yfDbWZ5195E5duazE_
#define DEFINED_TYPEDEF_FOR_struct_nKwr8yfDbWZ5195E5duazE_

typedef struct
{
    real_T SIMULATION_TIME;
    boolean_T SCOPES_ALL;
    boolean_T SCOPES_EXT_WRENCHES;
    boolean_T SCOPES_GAIN_SCHE_INFO;
    boolean_T SCOPES_MAIN;
    boolean_T SCOPES_QP;
    boolean_T CHECK_LIMITS;
    boolean_T SAVE_WORKSPACE;
    real_T Ts;
    boolean_T ON_GAZEBO;
    boolean_T USE_MOTOR_REFLECTED_INERTIA;
    boolean_T INCLUDE_COUPLING;
    boolean_T USE_IMU4EST_BASE;
    boolean_T FILTER_IMU_YAW;
    boolean_T FILTER_IMU_PITCH;
    boolean_T CORRECT_NECK_IMU;
    boolean_T USE_QP_SOLVER;
    real_T LEFT_RIGHT_FOOT_IN_CONTACT[2];
    boolean_T LEFT_FOOT_IN_CONTACT_AT_0;
    boolean_T DEMO_MOVEMENTS;
    boolean_T SMOOTH_COM_DES;
    boolean_T SMOOTH_JOINT_DES;
    real_T noOscillationTime;
    real_T directionOfOscillation[3];
    real_T amplitudeOfOscillation;
    real_T frequencyOfOscillation;
    real_T Gamma[529];
    real_T I_m[529];
    real_T T[529];
    real_T K_ff;
} struct_nKwr8yfDbWZ5195E5duazE;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_xTYJupe8umqLUY4vWPaUnH_
#define DEFINED_TYPEDEF_FOR_struct_xTYJupe8umqLUY4vWPaUnH_

typedef struct
{
    real_T SM_MASK_COORDINATOR;
    real_T SM_MASK_YOGA;
    real_T SM_TYPE_BIN;
    real_T smoothingTimeCoM_Joints[13];
    real_T joints_pauseBetweenYogaMoves;
    real_T wrench_thresholdContactOn;
    real_T wrench_thresholdContactOff;
    real_T CoM_threshold;
    real_T joints_thresholdNotInContact;
    real_T joints_thresholdInContact;
    real_T stateAt0;
    real_T CoM_delta[39];
    real_T tBalancing;
    real_T tBalancingBeforeYoga;
    boolean_T yogaExtended;
    boolean_T skipYoga;
    boolean_T demoOnlyBalancing;
    boolean_T demoStartsOnRightSupport;
    boolean_T yogaAlsoOnRightFoot;
    boolean_T yogaInLoop;
    real_T joints_references[299];
    real_T joints_leftYogaRef[624];
    real_T joints_rightYogaRef[624];
} struct_xTYJupe8umqLUY4vWPaUnH;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_iWQaDqPVqHoUwXRopcYWkG_
#define DEFINED_TYPEDEF_FOR_struct_iWQaDqPVqHoUwXRopcYWkG_

typedef struct
{
    real_T pinvDamp_nu_b;
    real_T pinvDamp;
    real_T pinvTol;
    real_T impedances;
    real_T dampings;
    real_T HessianQP;
} struct_iWQaDqPVqHoUwXRopcYWkG;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_r750FiCaocHqTuxb1nF3DC_
#define DEFINED_TYPEDEF_FOR_struct_r750FiCaocHqTuxb1nF3DC_

typedef struct
{
    real_T KP_COM[39];
    real_T KD_COM[39];
    real_T KP_AngularMomentum;
    real_T KD_AngularMomentum;
    real_T impedances[299];
    real_T dampings[23];
    real_T SmoothingTimeGainScheduling;
} struct_r750FiCaocHqTuxb1nF3DC;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_5rbQfDr6lxmxJDTJAmE7vH_
#define DEFINED_TYPEDEF_FOR_struct_5rbQfDr6lxmxJDTJAmE7vH_

typedef struct
{
    real_T torque;
} struct_5rbQfDr6lxmxJDTJAmE7vH;

#endif

/* Parameters (auto storage) */
typedef struct P_torqueBalancingYoga_T_ P_torqueBalancingYoga_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_torqueBalancingYoga_T RT_MODEL_torqueBalancingYoga_T;

#endif /* RTW_HEADER_torqueBalancingYoga_types_h_ */
