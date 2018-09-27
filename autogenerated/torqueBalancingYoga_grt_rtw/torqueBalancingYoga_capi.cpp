/*
 * torqueBalancingYoga_capi.cpp
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueBalancingYoga".
 *
 * Model version              : 1.3268
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C++ source code generated on : Fri Aug 10 20:57:30 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objective: Debugging
 * Validation result: Not run
 */

#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "torqueBalancingYoga_capi_host.h"
#define sizeof(s)((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s, el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)
#else /* HOST_CAPI_BUILD */
#include "builtin_typeid_types.h"
#include "torqueBalancingYoga.h"
#include "torqueBalancingYoga_capi.h"
#include "torqueBalancingYoga_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s) (NULL)
#else
#define TARGET_CONST const
#define TARGET_STRING(s) (s)
#endif
#endif /* HOST_CAPI_BUILD */

/* Block output signal information */
static rtwCAPI_Signals rtBlockSignals[] = {
    /* addrMapIndex, sysNum, blockPath,
     * signalName, portNumber, dataTypeIndex, dimIndex, fxpIndex, sTimeIndex
     */
    {0, 0, TARGET_STRING("torqueBalancingYoga/Saturation"), TARGET_STRING(""), 0, 0, 0, 0, 0},

    {1,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     0,
     0,
     1,
     0,
     0},

    {2,
     8,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     1,
     0,
     0,
     0,
     0},

    {3,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     2,
     0,
     2,
     0,
     1},

    {4,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     3,
     0,
     0,
     0,
     1},

    {5,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     4,
     0,
     3,
     0,
     1},

    {6,
     6,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING("nu_b"),
     5,
     0,
     4,
     0,
     0},

    {7,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     5,
     0,
     0,
     0,
     0},

    {8,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     6,
     0,
     0,
     0,
     0},

    {9,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     7,
     0,
     3,
     0,
     1},

    {10,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     8,
     0,
     5,
     0,
     1},

    {11,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     9,
     0,
     5,
     0,
     1},

    {12,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator"),
     TARGET_STRING(""),
     10,
     0,
     6,
     0,
     0},

    {13,
     24,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING("CoM_des"),
     0,
     0,
     5,
     0,
     0},

    {14,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     1},

    {15,
     24,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     1,
     0,
     0,
     0,
     0},

    {16,
     24,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     2,
     0,
     2,
     0,
     0},

    {17,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     3,
     0,
     7,
     0,
     0},

    {18,
     24,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     4,
     0,
     3,
     0,
     0},

    {19,
     12,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING("nu_b"),
     5,
     0,
     4,
     0,
     0},

    {20,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     5,
     0,
     0,
     0,
     0},

    {21,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     6,
     0,
     0,
     0,
     0},

    {22,
     24,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     7,
     0,
     3,
     0,
     0},

    {23,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     8,
     0,
     8,
     0,
     0},

    {24,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     9,
     0,
     8,
     0,
     0},

    {25,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga"),
     TARGET_STRING(""),
     10,
     0,
     6,
     0,
     0},

    {26,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Constant2"),
     TARGET_STRING(""),
     0,
     1,
     3,
     0,
     1},

    {27,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Multiport Switch1"),
     TARGET_STRING(""),
     0,
     0,
     9,
     0,
     0},

    {28,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     4,
     0,
     10,
     0,
     0},

    {29,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     5,
     0,
     4,
     0,
     0},

    {30,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     6,
     0,
     11,
     0,
     0},

    {31,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     7,
     0,
     12,
     0,
     0},

    {32,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     8,
     0,
     13,
     0,
     0},

    {33,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     9,
     0,
     14,
     0,
     0},

    {34,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     10,
     0,
     15,
     0,
     0},

    {35,
     26,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Balancing Controller "),
     TARGET_STRING(""),
     11,
     0,
     16,
     0,
     0},

    {36,
     45,
     TARGET_STRING("torqueBalancingYoga/emergency stop: joint limits/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     3,
     0,
     0},

    {37,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Compare To Zero/Compare"),
     TARGET_STRING(""),
     0,
     2,
     3,
     0,
     1},

    {38,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Compare To Zero1/Compare"),
     TARGET_STRING(""),
     0,
     2,
     3,
     0,
     1},

    {39,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Get Measurement/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {40,
     10,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Internal Coordinator/jointAngles"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {41,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant1"),
     TARGET_STRING(""),
     0,
     0,
     2,
     0,
     1},

    {42,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant3"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     1},

    {43,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant4"),
     TARGET_STRING(""),
     0,
     0,
     3,
     0,
     1},

    {44,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant5"),
     TARGET_STRING(""),
     0,
     0,
     5,
     0,
     1},

    {45,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant6"),
     TARGET_STRING(""),
     0,
     0,
     5,
     0,
     1},

    {46,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal "
                   "Coordinator/joints.smoothingTime"),
     TARGET_STRING(""),
     0,
     0,
     3,
     0,
     1},

    {47,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Reshape1"),
     TARGET_STRING(""),
     0,
     0,
     6,
     0,
     0},

    {48,
     10,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Internal Coordinator/IMU measurements"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {49,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Switch1"),
     TARGET_STRING(""),
     0,
     0,
     1,
     0,
     0},

    {50,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Select State and "
                   "References/Logical Operator1"),
     TARGET_STRING(""),
     0,
     2,
     3,
     0,
     1},

    {51,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Select State and "
                   "References/Minimum Jerk Trajectory Generator1"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {52,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Select State and "
                   "References/Minimum Jerk Trajectory Generator2"),
     TARGET_STRING(""),
     0,
     0,
     5,
     0,
     0},

    {53,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Select State and "
                   "References/Minimum Jerk Trajectory Generator2"),
     TARGET_STRING(""),
     1,
     0,
     5,
     0,
     0},

    {54,
     0,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Select State and "
                   "References/Minimum Jerk Trajectory Generator2"),
     TARGET_STRING(""),
     2,
     0,
     5,
     0,
     0},

    {55,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Select State and References/Switch5"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {56,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/jointAngles"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {57,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {58,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING(""),
     1,
     0,
     5,
     0,
     0},

    {59,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING(""),
     2,
     0,
     0,
     0,
     0},

    {60,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING("constraintsSmooth"),
     3,
     0,
     2,
     0,
     0},

    {61,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING(""),
     7,
     0,
     3,
     0,
     0},

    {62,
     24,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/stateMachineYogaFCN"),
     TARGET_STRING(""),
     8,
     0,
     3,
     0,
     0},

    {63,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Constant1"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     1},

    {64,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Reshape1"),
     TARGET_STRING(""),
     0,
     0,
     6,
     0,
     0},

    {65,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Minimum Jerk "
                   "Trajectory Generator"),
     TARGET_STRING(""),
     0,
     0,
     18,
     0,
     0},

    {66,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/inertial"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {67,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/left_foot_wrench"),
     TARGET_STRING("wL_WBDT"),
     0,
     0,
     4,
     0,
     0},

    {68,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/right_foot_wrench"),
     TARGET_STRING("wR_WBDT"),
     0,
     0,
     4,
     0,
     0},

    {69,
     37,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/References for H"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {70,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/inertial"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {71,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Sum"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {72,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Switch"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {73,
     44,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques (motor reflected "
                   "inertia)/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}"),
     TARGET_STRING(""),
     0,
     0,
     19,
     0,
     1},

    {74,
     46,
     TARGET_STRING("torqueBalancingYoga/emergency stop: joint limits/Get Limits/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {75,
     46,
     TARGET_STRING("torqueBalancingYoga/emergency stop: joint limits/Get Limits/S-Function"),
     TARGET_STRING(""),
     1,
     0,
     0,
     0,
     0},

    {76,
     50,
     TARGET_STRING("torqueBalancingYoga/tauDot Saturation/holder /MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {77,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Centroidal Momentum/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {78,
     0,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Get Bias Forces/Gain"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {79,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Get Bias Forces/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     18,
     0,
     0},

    {80,
     0,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Mass Matrix/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     20,
     0,
     0},

    {81,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Kinematics/DotJ Nu l_sole /S-Function"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {82,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Kinematics/DotJ Nu r_sole   /S-Function"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {83,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian com/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {84,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian l_sole/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {85,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/Dynamics and Kinematics/Kinematics/Jacobian r_sole/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {86,
     0,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Kinematics/l_sole/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {87,
     0,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Kinematics/r_sole/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {88,
     6,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute "
                   "State Velocity/Compute Base Velocity"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {89,
     8,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/holder "
                   "/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {90,
     9,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/holder "
                   "1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     5,
     0,
     0},

    {91,
     12,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute "
                   "State Velocity/Compute Base Velocity"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {92,
     22,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/holder "
                   "1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     5,
     0,
     0},

    {93,
     23,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/holder "
                   "2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {94,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Centroidal "
                   "Momentum/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {95,
     0,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Jacobian/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {96,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum "
                   "integral/Jacobian1/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {97,
     41,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {98,
     41,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot"),
     TARGET_STRING(""),
     1,
     0,
     3,
     0,
     0},

    {99,
     42,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet"),
     TARGET_STRING("QPBlock_TwoFeet"),
     0,
     0,
     14,
     0,
     0},

    {100,
     42,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet"),
     TARGET_STRING(""),
     1,
     0,
     3,
     0,
     0},

    {101,
     39,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/ContactsTransition"),
     TARGET_STRING(""),
     0,
     2,
     3,
     0,
     0},

    {102,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/not"),
     TARGET_STRING(""),
     0,
     2,
     3,
     0,
     0},

    {103,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques (motor reflected "
                   "inertia)/Get Measurement/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {104,
     0,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Kinematics/xCom/CoM/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {105,
     2,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /fromImuToHomogeousTransformFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {106,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /Neck Position"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {107,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /Switch6"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {108,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {109,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/RFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {110,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute "
                   "State Velocity/Get Measurement/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {111,
     10,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Internal Coordinator/xCom/CoM/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {112,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute "
                   "State Velocity/Get Measurement/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     0,
     0,
     0},

    {113,
     13,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /fromImuToHomogeousTransformFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {114,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to base link transform /Neck Position"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {115,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to base link transform /Switch6"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {116,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {117,
     17,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/RFoot to base link transform/fromImuToHomogeousTransformFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {118,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/Neck Position"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {119,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/Switch6"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {120,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {121,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom/CoM/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {122,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/xCom1/CoM/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {123,
     29,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /fromImuToHomogeousTransformFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {124,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/LFoot to base link transform /Neck Position"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {125,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/LFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {126,
     33,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/fromImuToHomogeousTransformFCN"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {127,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/RFoot to base link transform/Neck Position"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {128,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/RFoot to world transform (fixed base)/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {129,
     41,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot/QP"),
     TARGET_STRING(""),
     0,
     0,
     4,
     0,
     0},

    {130,
     41,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/One Foot/QP"),
     TARGET_STRING(""),
     1,
     0,
     3,
     0,
     0},

    {131,
     42,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet/QP"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {132,
     42,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute joint torques/QPSolver/Two Feet/QP"),
     TARGET_STRING(""),
     1,
     0,
     3,
     0,
     0},

    {133,
     10,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
         "fixed_link/LFoot to base link transform /Fixed base to imu transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {134,
     10,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
         "fixed_link/LFoot to base link transform /Fixed base to root link transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {135,
     3,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /holder 1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {136,
     4,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /holder 2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {137,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute "
                   "State Velocity/Feet Jacobians/Jacobian LFoot/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {138,
     10,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Compute "
                   "State Velocity/Feet Jacobians/Jacobian RFoot/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {139,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute "
                   "State Velocity/Feet Jacobians/Jacobian LFoot/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {140,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute "
                   "State Velocity/Feet Jacobians/Jacobian RFoot/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     21,
     0,
     0},

    {141,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /Fixed base to imu transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {142,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to base link transform /Fixed base to root link "
                   "transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {143,
     14,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /holder 1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {144,
     15,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /holder 2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {145,
     25,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/RFoot to base link transform/Fixed base to imu transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {146,
     25,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/Fixed base to root link "
                   "transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {147,
     18,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/holder 1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {148,
     19,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/holder 2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {149,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/LFoot to base link transform /Fixed base to imu "
                   "transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {150,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/LFoot to base link transform /Fixed base to root "
                   "link transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {151,
     30,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /holder 1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {152,
     31,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /holder 2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {153,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/RFoot to base link transform/Fixed base to imu "
                   "transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {154,
     0,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/RFoot to base link transform/Fixed base to root "
                   "link transform/S-Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {155,
     34,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/holder 1/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     17,
     0,
     0},

    {156,
     35,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/holder 2/MATLAB Function"),
     TARGET_STRING(""),
     0,
     0,
     14,
     0,
     0},

    {0, 0, (NULL), (NULL), 0, 0, 0, 0, 0}};

static rtwCAPI_BlockParameters rtBlockParameters[] = {
    /* addrMapIndex, blockPath,
     * paramName, dataTypeIndex, dimIndex, fixPtIdx
     */
    {157,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Constant"),
     TARGET_STRING("Value"),
     1,
     3,
     0},

    {158,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Constant2"),
     TARGET_STRING("Value"),
     1,
     3,
     0},

    {159,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Coordinator"),
     TARGET_STRING("BitMask"),
     1,
     3,
     0},

    {160,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Yoga"),
     TARGET_STRING("BitMask"),
     1,
     3,
     0},

    {161,
     TARGET_STRING("torqueBalancingYoga/emergency stop: joint limits/index1"),
     TARGET_STRING("Value"),
     0,
     3,
     0},

    {162,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Compare To Zero/Constant"),
     TARGET_STRING("Value"),
     1,
     3,
     0},

    {163,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Compare To Zero1/Constant"),
     TARGET_STRING("Value"),
     1,
     3,
     0},

    {164,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant2"),
     TARGET_STRING("Value"),
     0,
     4,
     0},

    {165,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant3"),
     TARGET_STRING("Value"),
     0,
     7,
     0},

    {166,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant4"),
     TARGET_STRING("Value"),
     0,
     3,
     0},

    {167,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant5"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {168,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Constant6"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {169,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Constant1"),
     TARGET_STRING("Value"),
     0,
     4,
     0},

    {170,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Switch"),
     TARGET_STRING("Threshold"),
     0,
     3,
     0},

    {171,
     TARGET_STRING("torqueBalancingYoga/tauDot Saturation/holder /Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {172,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Get Bias Forces/Constant"),
     TARGET_STRING("Value"),
     0,
     4,
     0},

    {173,
     TARGET_STRING("torqueBalancingYoga/Dynamics and Kinematics/Dynamics/Get Bias Forces/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {174,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/Constant7"),
     TARGET_STRING("Value"),
     0,
     17,
     0},

    {175,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/holder "
                   "/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {176,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/holder "
                   "1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {177,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/Constant7"),
     TARGET_STRING("Value"),
     0,
     17,
     0},

    {178,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/holder "
                   "1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {179,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/holder "
                   "2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {180,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/Constant7"),
     TARGET_STRING("Value"),
     0,
     17,
     0},

    {181,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /holder 1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {182,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /holder 2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {183,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /neck transform/Constant"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {184,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/Internal Coordinator/Base to "
                   "fixed_link/LFoot to base link transform /neck transform/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {185,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /holder 1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {186,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/LFoot to base link transform /holder 2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {187,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to base link transform /neck transform/Constant"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {188,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/LFoot to base link transform /neck transform/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {189,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/RFoot to base link transform/holder 1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {190,
     TARGET_STRING(
         "torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base to fixed "
         "link transform/RFoot to base link transform/holder 2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {191,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/neck transform/Constant"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {192,
     TARGET_STRING("torqueBalancingYoga/Robot State and References/State Machine Yoga/Compute base "
                   "to fixed link transform/RFoot to base link transform/neck transform/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {193,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /holder 1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {194,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /holder 2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {195,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /neck transform/Constant"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {196,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/LFoot to base link transform /neck transform/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {197,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/holder 1/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {198,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/holder 2/Compare To Constant"),
     TARGET_STRING("const"),
     0,
     3,
     0},

    {199,
     TARGET_STRING(
         "torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute base to "
         "fixed link transform/RFoot to base link transform/neck transform/Constant"),
     TARGET_STRING("Value"),
     0,
     5,
     0},

    {200,
     TARGET_STRING("torqueBalancingYoga/controller_QP/Compute angular momentum integral/Compute "
                   "base to fixed link transform/RFoot to base link transform/neck transform/Gain"),
     TARGET_STRING("Gain"),
     0,
     3,
     0},

    {0, (NULL), (NULL), 0, 0, 0}};

/* Tunable variable parameters */
static rtwCAPI_ModelParameters rtModelParameters[] = {
    /* addrMapIndex, varName, dataTypeIndex, dimIndex, fixPtIndex */
    {201, TARGET_STRING("Sm"), 3, 3, 0},

    {202, TARGET_STRING("Config"), 4, 3, 0},

    {203, TARGET_STRING("Gain"), 5, 3, 0},

    {204, TARGET_STRING("Reg"), 6, 3, 0},

    {205, TARGET_STRING("ConstraintsMatrix"), 0, 11, 0},

    {206, TARGET_STRING("ROBOT_DOF_FOR_SIMULINK"), 0, 19, 0},

    {207, TARGET_STRING("bVectorConstraints"), 0, 12, 0},

    {208, TARGET_STRING("Sat"), 7, 3, 0},

    {0, (NULL), 0, 0, 0}};

#ifndef HOST_CAPI_BUILD

/* Initialize Data Address */
static void torqueBalancingYoga_InitializeDataAddr(void* dataAddr[],
                                                   B_torqueBalancingYoga_T* torqueBalancingYoga_B,
                                                   P_torqueBalancingYoga_T* torqueBalancingYoga_P)
{
    dataAddr[0] = (void*) (&torqueBalancingYoga_B->Saturation[0]);
    dataAddr[1] = (void*) (&torqueBalancingYoga_B->Switch1[0]);
    dataAddr[2] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_k.s0[0]);
    dataAddr[3] = (void*) (&torqueBalancingYoga_B->Constant1_c[0]);
    dataAddr[4] = (void*) (&torqueBalancingYoga_B->Constant3[0]);
    dataAddr[5] = (void*) (&torqueBalancingYoga_B->Constant4);
    dataAddr[6] = (void*) (&torqueBalancingYoga_B->nu_b_i[0]);
    dataAddr[7] = (void*) (&torqueBalancingYoga_B->SFunction_m[0]);
    dataAddr[8] = (void*) (&torqueBalancingYoga_B->jointAngles_p[0]);
    dataAddr[9] = (void*) (&torqueBalancingYoga_B->jointssmoothingTime);
    dataAddr[10] = (void*) (&torqueBalancingYoga_B->Constant5[0]);
    dataAddr[11] = (void*) (&torqueBalancingYoga_B->Constant6[0]);
    dataAddr[12] = (void*) (&torqueBalancingYoga_B->Reshape1_p[0]);
    dataAddr[13] = (void*) (&torqueBalancingYoga_B->CoM_des[0]);
    dataAddr[14] = (void*) (&torqueBalancingYoga_B->Constant1[0]);
    dataAddr[15] = (void*) (&torqueBalancingYoga_B->qj_des[0]);
    dataAddr[16] = (void*) (&torqueBalancingYoga_B->constraints[0]);
    dataAddr[17] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator[0]);
    dataAddr[18] = (void*) (&torqueBalancingYoga_B->currentState);
    dataAddr[19] = (void*) (&torqueBalancingYoga_B->nu_b[0]);
    dataAddr[20] = (void*) (&torqueBalancingYoga_B->SFunction_nr[0]);
    dataAddr[21] = (void*) (&torqueBalancingYoga_B->jointAngles[0]);
    dataAddr[22] = (void*) (&torqueBalancingYoga_B->jointsSmoothingTime);
    dataAddr[23] = (void*) (((&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator[0]) + 23));
    dataAddr[24] = (void*) (((&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator[0]) + 26));
    dataAddr[25] = (void*) (&torqueBalancingYoga_B->Reshape1[0]);
    dataAddr[26] = (void*) (&torqueBalancingYoga_B->Constant2);
    dataAddr[27] = (void*) (&torqueBalancingYoga_B->MultiportSwitch1[0]);
    dataAddr[28] = (void*) (&torqueBalancingYoga_B->HessianMatrixQP1Foot[0]);
    dataAddr[29] = (void*) (&torqueBalancingYoga_B->gradientQP1Foot[0]);
    dataAddr[30] = (void*) (&torqueBalancingYoga_B->ConstraintsMatrixQP1Foot[0]);
    dataAddr[31] = (void*) (&torqueBalancingYoga_B->bVectorConstraintsQp1Foot[0]);
    dataAddr[32] = (void*) (&torqueBalancingYoga_B->HessianMatrixQP2Feet[0]);
    dataAddr[33] = (void*) (&torqueBalancingYoga_B->gradientQP2Feet[0]);
    dataAddr[34] = (void*) (&torqueBalancingYoga_B->ConstraintsMatrixQP2Feet[0]);
    dataAddr[35] = (void*) (&torqueBalancingYoga_B->bVectorConstraintsQp2Feet[0]);
    dataAddr[36] = (void*) (&torqueBalancingYoga_B->inRange);
    dataAddr[37] = (void*) (&torqueBalancingYoga_B->Compare_c);
    dataAddr[38] = (void*) (&torqueBalancingYoga_B->Compare);
    dataAddr[39] = (void*) (&torqueBalancingYoga_B->SFunction_d[0]);
    dataAddr[40] = (void*) (&torqueBalancingYoga_B->jointAngles_p[0]);
    dataAddr[41] = (void*) (&torqueBalancingYoga_B->Constant1_c[0]);
    dataAddr[42] = (void*) (&torqueBalancingYoga_B->Constant3[0]);
    dataAddr[43] = (void*) (&torqueBalancingYoga_B->Constant4);
    dataAddr[44] = (void*) (&torqueBalancingYoga_B->Constant5[0]);
    dataAddr[45] = (void*) (&torqueBalancingYoga_B->Constant6[0]);
    dataAddr[46] = (void*) (&torqueBalancingYoga_B->jointssmoothingTime);
    dataAddr[47] = (void*) (&torqueBalancingYoga_B->Reshape1_p[0]);
    dataAddr[48] = (void*) (&torqueBalancingYoga_B->IMUmeasurements[0]);
    dataAddr[49] = (void*) (&torqueBalancingYoga_B->Switch1[0]);
    dataAddr[50] = (void*) (&torqueBalancingYoga_B->LogicalOperator1);
    dataAddr[51] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator1[0]);
    dataAddr[52] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator2[0]);
    dataAddr[53] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerato_d[0]);
    dataAddr[54] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerato_b[0]);
    dataAddr[55] = (void*) (&torqueBalancingYoga_B->Switch5[0]);
    dataAddr[56] = (void*) (&torqueBalancingYoga_B->jointAngles[0]);
    dataAddr[57] = (void*) (&torqueBalancingYoga_B->w_H_b[0]);
    dataAddr[58] = (void*) (&torqueBalancingYoga_B->CoM_des[0]);
    dataAddr[59] = (void*) (&torqueBalancingYoga_B->qj_des[0]);
    dataAddr[60] = (void*) (&torqueBalancingYoga_B->constraints[0]);
    dataAddr[61] = (void*) (&torqueBalancingYoga_B->currentState);
    dataAddr[62] = (void*) (&torqueBalancingYoga_B->jointsSmoothingTime);
    dataAddr[63] = (void*) (&torqueBalancingYoga_B->Constant1[0]);
    dataAddr[64] = (void*) (&torqueBalancingYoga_B->Reshape1[0]);
    dataAddr[65] = (void*) (&torqueBalancingYoga_B->MinimumJerkTrajectoryGenerator[0]);
    dataAddr[66] = (void*) (&torqueBalancingYoga_B->inertial_n[0]);
    dataAddr[67] = (void*) (&torqueBalancingYoga_B->wL_WBDT[0]);
    dataAddr[68] = (void*) (&torqueBalancingYoga_B->wR_WBDT[0]);
    dataAddr[69] = (void*) (&torqueBalancingYoga_B->nu_b_equivalent[0]);
    dataAddr[70] = (void*) (&torqueBalancingYoga_B->inertial[0]);
    dataAddr[71] = (void*) (&torqueBalancingYoga_B->Sum[0]);
    dataAddr[72] = (void*) (&torqueBalancingYoga_B->Switch[0]);
    dataAddr[73] = (void*) (&torqueBalancingYoga_B->reflectedInertia[0]);
    dataAddr[74] = (void*) (&torqueBalancingYoga_B->SFunction_o1[0]);
    dataAddr[75] = (void*) (&torqueBalancingYoga_B->SFunction_o2[0]);
    dataAddr[76] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_p.s0[0]);
    dataAddr[77] = (void*) (&torqueBalancingYoga_B->SFunction_c[0]);
    dataAddr[78] = (void*) (&torqueBalancingYoga_B->Gain[0]);
    dataAddr[79] = (void*) (&torqueBalancingYoga_B->SFunction_f[0]);
    dataAddr[80] = (void*) (&torqueBalancingYoga_B->SFunction_b[0]);
    dataAddr[81] = (void*) (&torqueBalancingYoga_B->SFunction_ej[0]);
    dataAddr[82] = (void*) (&torqueBalancingYoga_B->SFunction_h[0]);
    dataAddr[83] = (void*) (&torqueBalancingYoga_B->SFunction_n[0]);
    dataAddr[84] = (void*) (&torqueBalancingYoga_B->SFunction_c4[0]);
    dataAddr[85] = (void*) (&torqueBalancingYoga_B->SFunction_b3[0]);
    dataAddr[86] = (void*) (&torqueBalancingYoga_B->SFunction_oe[0]);
    dataAddr[87] = (void*) (&torqueBalancingYoga_B->SFunction_j[0]);
    dataAddr[88] = (void*) (&torqueBalancingYoga_B->nu_b_i[0]);
    dataAddr[89] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_k.s0[0]);
    dataAddr[90] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_a.s0[0]);
    dataAddr[91] = (void*) (&torqueBalancingYoga_B->nu_b[0]);
    dataAddr[92] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_f.s0[0]);
    dataAddr[93] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_d4.s0[0]);
    dataAddr[94] = (void*) (&torqueBalancingYoga_B->SFunction_on[0]);
    dataAddr[95] = (void*) (&torqueBalancingYoga_B->SFunction_e[0]);
    dataAddr[96] = (void*) (&torqueBalancingYoga_B->SFunction_cf[0]);
    dataAddr[97] = (void*) (&torqueBalancingYoga_B->QP_o1_f[0]);
    dataAddr[98] = (void*) (&torqueBalancingYoga_B->QP_o2_f);
    dataAddr[99] = (void*) (&torqueBalancingYoga_B->QP_o1[0]);
    dataAddr[100] = (void*) (&torqueBalancingYoga_B->QP_o2);
    dataAddr[101] = (void*) (&torqueBalancingYoga_B->onOneFoot);
    dataAddr[102] = (void*) (&torqueBalancingYoga_B->not_p);
    dataAddr[103] = (void*) (&torqueBalancingYoga_B->SFunction[0]);
    dataAddr[104] = (void*) (&torqueBalancingYoga_B->SFunction_dk[0]);
    dataAddr[105] = (void*) (&torqueBalancingYoga_B->sf_fromImuToHomogeousTransfor_f.w_H_b[0]);
    dataAddr[106] = (void*) (&torqueBalancingYoga_B->NeckPosition_c[0]);
    dataAddr[107] = (void*) (&torqueBalancingYoga_B->Switch6_b[0]);
    dataAddr[108] = (void*) (&torqueBalancingYoga_B->SFunction_k0[0]);
    dataAddr[109] = (void*) (&torqueBalancingYoga_B->SFunction_bdr[0]);
    dataAddr[110] = (void*) (&torqueBalancingYoga_B->SFunction_m[0]);
    dataAddr[111] = (void*) (&torqueBalancingYoga_B->SFunction_ak[0]);
    dataAddr[112] = (void*) (&torqueBalancingYoga_B->SFunction_nr[0]);
    dataAddr[113] = (void*) (&torqueBalancingYoga_B->sf_fromImuToHomogeousTransfo_fb.w_H_b[0]);
    dataAddr[114] = (void*) (&torqueBalancingYoga_B->NeckPosition_p[0]);
    dataAddr[115] = (void*) (&torqueBalancingYoga_B->Switch6[0]);
    dataAddr[116] = (void*) (&torqueBalancingYoga_B->SFunction_ao[0]);
    dataAddr[117] = (void*) (&torqueBalancingYoga_B->sf_fromImuToHomogeousTransfor_g.w_H_b[0]);
    dataAddr[118] = (void*) (&torqueBalancingYoga_B->NeckPosition_k[0]);
    dataAddr[119] = (void*) (&torqueBalancingYoga_B->Switch6_g[0]);
    dataAddr[120] = (void*) (&torqueBalancingYoga_B->SFunction_i[0]);
    dataAddr[121] = (void*) (&torqueBalancingYoga_B->SFunction_jj[0]);
    dataAddr[122] = (void*) (&torqueBalancingYoga_B->SFunction_dn[0]);
    dataAddr[123] = (void*) (&torqueBalancingYoga_B->sf_fromImuToHomogeousTransformF.w_H_b[0]);
    dataAddr[124] = (void*) (&torqueBalancingYoga_B->NeckPosition[0]);
    dataAddr[125] = (void*) (&torqueBalancingYoga_B->SFunction_br[0]);
    dataAddr[126] = (void*) (&torqueBalancingYoga_B->sf_fromImuToHomogeousTransfor_e.w_H_b[0]);
    dataAddr[127] = (void*) (&torqueBalancingYoga_B->NeckPosition_m[0]);
    dataAddr[128] = (void*) (&torqueBalancingYoga_B->SFunction_a[0]);
    dataAddr[129] = (void*) (&torqueBalancingYoga_B->QP_o1_f[0]);
    dataAddr[130] = (void*) (&torqueBalancingYoga_B->QP_o2_f);
    dataAddr[131] = (void*) (&torqueBalancingYoga_B->QP_o1[0]);
    dataAddr[132] = (void*) (&torqueBalancingYoga_B->QP_o2);
    dataAddr[133] = (void*) (&torqueBalancingYoga_B->SFunction_jo[0]);
    dataAddr[134] = (void*) (&torqueBalancingYoga_B->SFunction_l[0]);
    dataAddr[135] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_os.s0[0]);
    dataAddr[136] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_l5.s0[0]);
    dataAddr[137] = (void*) (&torqueBalancingYoga_B->SFunction_p[0]);
    dataAddr[138] = (void*) (&torqueBalancingYoga_B->SFunction_i0[0]);
    dataAddr[139] = (void*) (&torqueBalancingYoga_B->SFunction_e0[0]);
    dataAddr[140] = (void*) (&torqueBalancingYoga_B->SFunction_da[0]);
    dataAddr[141] = (void*) (&torqueBalancingYoga_B->SFunction_k[0]);
    dataAddr[142] = (void*) (&torqueBalancingYoga_B->SFunction_nj[0]);
    dataAddr[143] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_i.s0[0]);
    dataAddr[144] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_oa.s0[0]);
    dataAddr[145] = (void*) (&torqueBalancingYoga_B->SFunction_oz[0]);
    dataAddr[146] = (void*) (&torqueBalancingYoga_B->SFunction_ef[0]);
    dataAddr[147] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_ad.s0[0]);
    dataAddr[148] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_j.s0[0]);
    dataAddr[149] = (void*) (&torqueBalancingYoga_B->SFunction_bd[0]);
    dataAddr[150] = (void*) (&torqueBalancingYoga_B->SFunction_d2[0]);
    dataAddr[151] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction.s0[0]);
    dataAddr[152] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_h.s0[0]);
    dataAddr[153] = (void*) (&torqueBalancingYoga_B->SFunction_o[0]);
    dataAddr[154] = (void*) (&torqueBalancingYoga_B->SFunction_c2[0]);
    dataAddr[155] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_o.s0[0]);
    dataAddr[156] = (void*) (&torqueBalancingYoga_B->sf_MATLABFunction_l.s0[0]);
    dataAddr[157] = (void*) (&torqueBalancingYoga_P->Constant_Value_it);
    dataAddr[158] = (void*) (&torqueBalancingYoga_P->Constant2_Value_o);
    dataAddr[159] = (void*) (&torqueBalancingYoga_P->Coordinator_BitMask);
    dataAddr[160] = (void*) (&torqueBalancingYoga_P->Yoga_BitMask);
    dataAddr[161] = (void*) (&torqueBalancingYoga_P->index1_Value);
    dataAddr[162] = (void*) (&torqueBalancingYoga_P->Constant_Value_dx);
    dataAddr[163] = (void*) (&torqueBalancingYoga_P->Constant_Value_n);
    dataAddr[164] = (void*) (&torqueBalancingYoga_P->Constant2_Value[0]);
    dataAddr[165] = (void*) (&torqueBalancingYoga_P->Constant3_Value[0]);
    dataAddr[166] = (void*) (&torqueBalancingYoga_P->Constant4_Value);
    dataAddr[167] = (void*) (&torqueBalancingYoga_P->Constant5_Value[0]);
    dataAddr[168] = (void*) (&torqueBalancingYoga_P->Constant6_Value[0]);
    dataAddr[169] = (void*) (&torqueBalancingYoga_P->Constant1_Value[0]);
    dataAddr[170] = (void*) (&torqueBalancingYoga_P->Switch_Threshold);
    dataAddr[171] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_l);
    dataAddr[172] = (void*) (&torqueBalancingYoga_P->Constant_Value_d[0]);
    dataAddr[173] = (void*) (&torqueBalancingYoga_P->Gain_Gain_o);
    dataAddr[174] = (void*) (&torqueBalancingYoga_P->Constant7_Value[0]);
    dataAddr[175] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_ok);
    dataAddr[176] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_o);
    dataAddr[177] = (void*) (&torqueBalancingYoga_P->Constant7_Value_b[0]);
    dataAddr[178] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_ph);
    dataAddr[179] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_pf);
    dataAddr[180] = (void*) (&torqueBalancingYoga_P->Constant7_Value_o[0]);
    dataAddr[181] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const);
    dataAddr[182] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_j);
    dataAddr[183] = (void*) (&torqueBalancingYoga_P->Constant_Value[0]);
    dataAddr[184] = (void*) (&torqueBalancingYoga_P->Gain_Gain);
    dataAddr[185] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_p);
    dataAddr[186] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_od);
    dataAddr[187] = (void*) (&torqueBalancingYoga_P->Constant_Value_c[0]);
    dataAddr[188] = (void*) (&torqueBalancingYoga_P->Gain_Gain_a);
    dataAddr[189] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_f);
    dataAddr[190] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_e);
    dataAddr[191] = (void*) (&torqueBalancingYoga_P->Constant_Value_cj[0]);
    dataAddr[192] = (void*) (&torqueBalancingYoga_P->Gain_Gain_h);
    dataAddr[193] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_h);
    dataAddr[194] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_c);
    dataAddr[195] = (void*) (&torqueBalancingYoga_P->Constant_Value_o[0]);
    dataAddr[196] = (void*) (&torqueBalancingYoga_P->Gain_Gain_hk);
    dataAddr[197] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_ea);
    dataAddr[198] = (void*) (&torqueBalancingYoga_P->CompareToConstant_const_e5);
    dataAddr[199] = (void*) (&torqueBalancingYoga_P->Constant_Value_i[0]);
    dataAddr[200] = (void*) (&torqueBalancingYoga_P->Gain_Gain_c);
    dataAddr[201] = (void*) (&torqueBalancingYoga_P->Sm);
    dataAddr[202] = (void*) (&torqueBalancingYoga_P->Config);
    dataAddr[203] = (void*) (&torqueBalancingYoga_P->Gain);
    dataAddr[204] = (void*) (&torqueBalancingYoga_P->Reg);
    dataAddr[205] = (void*) (&torqueBalancingYoga_P->ConstraintsMatrix[0]);
    dataAddr[206] = (void*) (&torqueBalancingYoga_P->ROBOT_DOF_FOR_SIMULINK[0]);
    dataAddr[207] = (void*) (&torqueBalancingYoga_P->bVectorConstraints[0]);
    dataAddr[208] = (void*) (&torqueBalancingYoga_P->Sat);
}

#endif

/* Initialize Data Run-Time Dimension Buffer Address */
#ifndef HOST_CAPI_BUILD

static void torqueBalancingYoga_InitializeVarDimsAddr(int32_T* vardimsAddr[])
{
    vardimsAddr[0] = (NULL);
}

#endif

#ifndef HOST_CAPI_BUILD

/* Initialize logging function pointers */
static void torqueBalancingYoga_InitializeLoggingFunctions(RTWLoggingFcnPtr loggingPtrs[])
{
    loggingPtrs[0] = (NULL);
    loggingPtrs[1] = (NULL);
    loggingPtrs[2] = (NULL);
    loggingPtrs[3] = (NULL);
    loggingPtrs[4] = (NULL);
    loggingPtrs[5] = (NULL);
    loggingPtrs[6] = (NULL);
    loggingPtrs[7] = (NULL);
    loggingPtrs[8] = (NULL);
    loggingPtrs[9] = (NULL);
    loggingPtrs[10] = (NULL);
    loggingPtrs[11] = (NULL);
    loggingPtrs[12] = (NULL);
    loggingPtrs[13] = (NULL);
    loggingPtrs[14] = (NULL);
    loggingPtrs[15] = (NULL);
    loggingPtrs[16] = (NULL);
    loggingPtrs[17] = (NULL);
    loggingPtrs[18] = (NULL);
    loggingPtrs[19] = (NULL);
    loggingPtrs[20] = (NULL);
    loggingPtrs[21] = (NULL);
    loggingPtrs[22] = (NULL);
    loggingPtrs[23] = (NULL);
    loggingPtrs[24] = (NULL);
    loggingPtrs[25] = (NULL);
    loggingPtrs[26] = (NULL);
    loggingPtrs[27] = (NULL);
    loggingPtrs[28] = (NULL);
    loggingPtrs[29] = (NULL);
    loggingPtrs[30] = (NULL);
    loggingPtrs[31] = (NULL);
    loggingPtrs[32] = (NULL);
    loggingPtrs[33] = (NULL);
    loggingPtrs[34] = (NULL);
    loggingPtrs[35] = (NULL);
    loggingPtrs[36] = (NULL);
    loggingPtrs[37] = (NULL);
    loggingPtrs[38] = (NULL);
    loggingPtrs[39] = (NULL);
    loggingPtrs[40] = (NULL);
    loggingPtrs[41] = (NULL);
    loggingPtrs[42] = (NULL);
    loggingPtrs[43] = (NULL);
    loggingPtrs[44] = (NULL);
    loggingPtrs[45] = (NULL);
    loggingPtrs[46] = (NULL);
    loggingPtrs[47] = (NULL);
    loggingPtrs[48] = (NULL);
    loggingPtrs[49] = (NULL);
    loggingPtrs[50] = (NULL);
    loggingPtrs[51] = (NULL);
    loggingPtrs[52] = (NULL);
    loggingPtrs[53] = (NULL);
    loggingPtrs[54] = (NULL);
    loggingPtrs[55] = (NULL);
    loggingPtrs[56] = (NULL);
    loggingPtrs[57] = (NULL);
    loggingPtrs[58] = (NULL);
    loggingPtrs[59] = (NULL);
    loggingPtrs[60] = (NULL);
    loggingPtrs[61] = (NULL);
    loggingPtrs[62] = (NULL);
    loggingPtrs[63] = (NULL);
    loggingPtrs[64] = (NULL);
    loggingPtrs[65] = (NULL);
    loggingPtrs[66] = (NULL);
    loggingPtrs[67] = (NULL);
    loggingPtrs[68] = (NULL);
    loggingPtrs[69] = (NULL);
    loggingPtrs[70] = (NULL);
    loggingPtrs[71] = (NULL);
    loggingPtrs[72] = (NULL);
    loggingPtrs[73] = (NULL);
    loggingPtrs[74] = (NULL);
    loggingPtrs[75] = (NULL);
    loggingPtrs[76] = (NULL);
    loggingPtrs[77] = (NULL);
    loggingPtrs[78] = (NULL);
    loggingPtrs[79] = (NULL);
    loggingPtrs[80] = (NULL);
    loggingPtrs[81] = (NULL);
    loggingPtrs[82] = (NULL);
    loggingPtrs[83] = (NULL);
    loggingPtrs[84] = (NULL);
    loggingPtrs[85] = (NULL);
    loggingPtrs[86] = (NULL);
    loggingPtrs[87] = (NULL);
    loggingPtrs[88] = (NULL);
    loggingPtrs[89] = (NULL);
    loggingPtrs[90] = (NULL);
    loggingPtrs[91] = (NULL);
    loggingPtrs[92] = (NULL);
    loggingPtrs[93] = (NULL);
    loggingPtrs[94] = (NULL);
    loggingPtrs[95] = (NULL);
    loggingPtrs[96] = (NULL);
    loggingPtrs[97] = (NULL);
    loggingPtrs[98] = (NULL);
    loggingPtrs[99] = (NULL);
    loggingPtrs[100] = (NULL);
    loggingPtrs[101] = (NULL);
    loggingPtrs[102] = (NULL);
    loggingPtrs[103] = (NULL);
    loggingPtrs[104] = (NULL);
    loggingPtrs[105] = (NULL);
    loggingPtrs[106] = (NULL);
    loggingPtrs[107] = (NULL);
    loggingPtrs[108] = (NULL);
    loggingPtrs[109] = (NULL);
    loggingPtrs[110] = (NULL);
    loggingPtrs[111] = (NULL);
    loggingPtrs[112] = (NULL);
    loggingPtrs[113] = (NULL);
    loggingPtrs[114] = (NULL);
    loggingPtrs[115] = (NULL);
    loggingPtrs[116] = (NULL);
    loggingPtrs[117] = (NULL);
    loggingPtrs[118] = (NULL);
    loggingPtrs[119] = (NULL);
    loggingPtrs[120] = (NULL);
    loggingPtrs[121] = (NULL);
    loggingPtrs[122] = (NULL);
    loggingPtrs[123] = (NULL);
    loggingPtrs[124] = (NULL);
    loggingPtrs[125] = (NULL);
    loggingPtrs[126] = (NULL);
    loggingPtrs[127] = (NULL);
    loggingPtrs[128] = (NULL);
    loggingPtrs[129] = (NULL);
    loggingPtrs[130] = (NULL);
    loggingPtrs[131] = (NULL);
    loggingPtrs[132] = (NULL);
    loggingPtrs[133] = (NULL);
    loggingPtrs[134] = (NULL);
    loggingPtrs[135] = (NULL);
    loggingPtrs[136] = (NULL);
    loggingPtrs[137] = (NULL);
    loggingPtrs[138] = (NULL);
    loggingPtrs[139] = (NULL);
    loggingPtrs[140] = (NULL);
    loggingPtrs[141] = (NULL);
    loggingPtrs[142] = (NULL);
    loggingPtrs[143] = (NULL);
    loggingPtrs[144] = (NULL);
    loggingPtrs[145] = (NULL);
    loggingPtrs[146] = (NULL);
    loggingPtrs[147] = (NULL);
    loggingPtrs[148] = (NULL);
    loggingPtrs[149] = (NULL);
    loggingPtrs[150] = (NULL);
    loggingPtrs[151] = (NULL);
    loggingPtrs[152] = (NULL);
    loggingPtrs[153] = (NULL);
    loggingPtrs[154] = (NULL);
    loggingPtrs[155] = (NULL);
    loggingPtrs[156] = (NULL);
    loggingPtrs[157] = (NULL);
    loggingPtrs[158] = (NULL);
    loggingPtrs[159] = (NULL);
    loggingPtrs[160] = (NULL);
    loggingPtrs[161] = (NULL);
    loggingPtrs[162] = (NULL);
    loggingPtrs[163] = (NULL);
    loggingPtrs[164] = (NULL);
    loggingPtrs[165] = (NULL);
    loggingPtrs[166] = (NULL);
    loggingPtrs[167] = (NULL);
    loggingPtrs[168] = (NULL);
    loggingPtrs[169] = (NULL);
    loggingPtrs[170] = (NULL);
    loggingPtrs[171] = (NULL);
    loggingPtrs[172] = (NULL);
    loggingPtrs[173] = (NULL);
    loggingPtrs[174] = (NULL);
    loggingPtrs[175] = (NULL);
    loggingPtrs[176] = (NULL);
    loggingPtrs[177] = (NULL);
    loggingPtrs[178] = (NULL);
    loggingPtrs[179] = (NULL);
    loggingPtrs[180] = (NULL);
    loggingPtrs[181] = (NULL);
    loggingPtrs[182] = (NULL);
    loggingPtrs[183] = (NULL);
    loggingPtrs[184] = (NULL);
    loggingPtrs[185] = (NULL);
    loggingPtrs[186] = (NULL);
    loggingPtrs[187] = (NULL);
    loggingPtrs[188] = (NULL);
    loggingPtrs[189] = (NULL);
    loggingPtrs[190] = (NULL);
    loggingPtrs[191] = (NULL);
    loggingPtrs[192] = (NULL);
    loggingPtrs[193] = (NULL);
    loggingPtrs[194] = (NULL);
    loggingPtrs[195] = (NULL);
    loggingPtrs[196] = (NULL);
    loggingPtrs[197] = (NULL);
    loggingPtrs[198] = (NULL);
    loggingPtrs[199] = (NULL);
    loggingPtrs[200] = (NULL);
    loggingPtrs[201] = (NULL);
    loggingPtrs[202] = (NULL);
    loggingPtrs[203] = (NULL);
    loggingPtrs[204] = (NULL);
    loggingPtrs[205] = (NULL);
    loggingPtrs[206] = (NULL);
    loggingPtrs[207] = (NULL);
    loggingPtrs[208] = (NULL);
}

#endif

/* Data Type Map - use dataTypeMapIndex to access this structure */
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap[] = {
    /* cName, mwName, numElements, elemMapIndex, dataSize, slDataId, *
     * isComplex, isPointer */
    {"double", "real_T", 0, 0, sizeof(real_T), SS_DOUBLE, 0, 0},

    {"unsigned char", "uint8_T", 0, 0, sizeof(uint8_T), SS_UINT8, 0, 0},

    {"unsigned char", "boolean_T", 0, 0, sizeof(boolean_T), SS_BOOLEAN, 0, 0},

    {"struct",
     "struct_xTYJupe8umqLUY4vWPaUnH",
     23,
     1,
     sizeof(struct_xTYJupe8umqLUY4vWPaUnH),
     SS_STRUCT,
     0,
     0},

    {"struct",
     "struct_hu1Fw0wiIZUzALhhiDLesH",
     32,
     24,
     sizeof(struct_hu1Fw0wiIZUzALhhiDLesH),
     SS_STRUCT,
     0,
     0},

    {"struct",
     "struct_r750FiCaocHqTuxb1nF3DC",
     7,
     56,
     sizeof(struct_r750FiCaocHqTuxb1nF3DC),
     SS_STRUCT,
     0,
     0},

    {"struct",
     "struct_iWQaDqPVqHoUwXRopcYWkG",
     6,
     63,
     sizeof(struct_iWQaDqPVqHoUwXRopcYWkG),
     SS_STRUCT,
     0,
     0},

    {"struct",
     "struct_5rbQfDr6lxmxJDTJAmE7vH",
     1,
     69,
     sizeof(struct_5rbQfDr6lxmxJDTJAmE7vH),
     SS_STRUCT,
     0,
     0}};

#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif

/* Structure Element Map - use elemMapIndex to access this structure */
static TARGET_CONST rtwCAPI_ElementMap rtElementMap[] = {
    /* elementName, elementOffset, dataTypeIndex, dimIndex, fxpIndex */
    {(NULL), 0, 0, 0, 0},

    {"SM_MASK_COORDINATOR",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, SM_MASK_COORDINATOR),
     0,
     22,
     0},

    {"SM_MASK_YOGA", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, SM_MASK_YOGA), 0, 22, 0},

    {"SM_TYPE_BIN", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, SM_TYPE_BIN), 0, 22, 0},

    {"smoothingTimeCoM_Joints",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, smoothingTimeCoM_Joints),
     0,
     23,
     0},

    {"joints_pauseBetweenYogaMoves",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_pauseBetweenYogaMoves),
     0,
     22,
     0},

    {"wrench_thresholdContactOn",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, wrench_thresholdContactOn),
     0,
     22,
     0},

    {"wrench_thresholdContactOff",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, wrench_thresholdContactOff),
     0,
     22,
     0},

    {"CoM_threshold", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, CoM_threshold), 0, 22, 0},

    {"joints_thresholdNotInContact",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_thresholdNotInContact),
     0,
     22,
     0},

    {"joints_thresholdInContact",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_thresholdInContact),
     0,
     22,
     0},

    {"stateAt0", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, stateAt0), 0, 22, 0},

    {"CoM_delta", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, CoM_delta), 0, 24, 0},

    {"tBalancing", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, tBalancing), 0, 22, 0},

    {"tBalancingBeforeYoga",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, tBalancingBeforeYoga),
     0,
     22,
     0},

    {"yogaExtended", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, yogaExtended), 2, 22, 0},

    {"skipYoga", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, skipYoga), 2, 22, 0},

    {"demoOnlyBalancing", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, demoOnlyBalancing), 2, 22, 0},

    {"demoStartsOnRightSupport",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, demoStartsOnRightSupport),
     2,
     22,
     0},

    {"yogaAlsoOnRightFoot",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, yogaAlsoOnRightFoot),
     2,
     22,
     0},

    {"yogaInLoop", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, yogaInLoop), 2, 22, 0},

    {"joints_references", rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_references), 0, 25, 0},

    {"joints_leftYogaRef",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_leftYogaRef),
     0,
     26,
     0},

    {"joints_rightYogaRef",
     rt_offsetof(struct_xTYJupe8umqLUY4vWPaUnH, joints_rightYogaRef),
     0,
     26,
     0},

    {"SIMULATION_TIME", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SIMULATION_TIME), 0, 22, 0},

    {"SCOPES_ALL", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SCOPES_ALL), 2, 22, 0},

    {"SCOPES_EXT_WRENCHES",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SCOPES_EXT_WRENCHES),
     2,
     22,
     0},

    {"SCOPES_GAIN_SCHE_INFO",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SCOPES_GAIN_SCHE_INFO),
     2,
     22,
     0},

    {"SCOPES_MAIN", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SCOPES_MAIN), 2, 22, 0},

    {"SCOPES_QP", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SCOPES_QP), 2, 22, 0},

    {"CHECK_LIMITS", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, CHECK_LIMITS), 2, 22, 0},

    {"SAVE_WORKSPACE", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SAVE_WORKSPACE), 2, 22, 0},

    {"Ts", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, Ts), 0, 22, 0},

    {"ON_GAZEBO", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, ON_GAZEBO), 2, 22, 0},

    {"SATURATE_TORQUE_DERIVATIVE",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SATURATE_TORQUE_DERIVATIVE),
     2,
     22,
     0},

    {"USE_MOTOR_REFLECTED_INERTIA",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, USE_MOTOR_REFLECTED_INERTIA),
     2,
     22,
     0},

    {"INCLUDE_COUPLING", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, INCLUDE_COUPLING), 2, 22, 0},

    {"USE_IMU4EST_BASE", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, USE_IMU4EST_BASE), 2, 22, 0},

    {"FILTER_IMU_YAW", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, FILTER_IMU_YAW), 2, 22, 0},

    {"FILTER_IMU_PITCH", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, FILTER_IMU_PITCH), 2, 22, 0},

    {"CORRECT_NECK_IMU", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, CORRECT_NECK_IMU), 2, 22, 0},

    {"USE_QP_SOLVER", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, USE_QP_SOLVER), 2, 22, 0},

    {"LEFT_RIGHT_FOOT_IN_CONTACT",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, LEFT_RIGHT_FOOT_IN_CONTACT),
     0,
     27,
     0},

    {"LEFT_FOOT_IN_CONTACT_AT_0",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, LEFT_FOOT_IN_CONTACT_AT_0),
     2,
     22,
     0},

    {"DEMO_MOVEMENTS", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, DEMO_MOVEMENTS), 2, 22, 0},

    {"SMOOTH_COM_DES", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SMOOTH_COM_DES), 2, 22, 0},

    {"SMOOTH_JOINT_DES", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, SMOOTH_JOINT_DES), 2, 22, 0},

    {"tauDot_maxAbs", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, tauDot_maxAbs), 0, 22, 0},

    {"noOscillationTime", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, noOscillationTime), 0, 22, 0},

    {"directionOfOscillation",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, directionOfOscillation),
     0,
     28,
     0},

    {"amplitudeOfOscillation",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, amplitudeOfOscillation),
     0,
     22,
     0},

    {"frequencyOfOscillation",
     rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, frequencyOfOscillation),
     0,
     22,
     0},

    {"Gamma", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, Gamma), 0, 19, 0},

    {"I_m", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, I_m), 0, 19, 0},

    {"T", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, T), 0, 19, 0},

    {"K_ff", rt_offsetof(struct_hu1Fw0wiIZUzALhhiDLesH, K_ff), 0, 22, 0},

    {"KP_COM", rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, KP_COM), 0, 24, 0},

    {"KD_COM", rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, KD_COM), 0, 24, 0},

    {"KP_AngularMomentum",
     rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, KP_AngularMomentum),
     0,
     22,
     0},

    {"KD_AngularMomentum",
     rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, KD_AngularMomentum),
     0,
     22,
     0},

    {"impedances", rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, impedances), 0, 25, 0},

    {"dampings", rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, dampings), 0, 29, 0},

    {"SmoothingTimeGainScheduling",
     rt_offsetof(struct_r750FiCaocHqTuxb1nF3DC, SmoothingTimeGainScheduling),
     0,
     22,
     0},

    {"pinvDamp_nu_b", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, pinvDamp_nu_b), 0, 22, 0},

    {"pinvDamp", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, pinvDamp), 0, 22, 0},

    {"pinvTol", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, pinvTol), 0, 22, 0},

    {"impedances", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, impedances), 0, 22, 0},

    {"dampings", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, dampings), 0, 22, 0},

    {"HessianQP", rt_offsetof(struct_iWQaDqPVqHoUwXRopcYWkG, HessianQP), 0, 22, 0},

    {"torque", rt_offsetof(struct_5rbQfDr6lxmxJDTJAmE7vH, torque), 0, 22, 0}};

/* Dimension Map - use dimensionMapIndex to access elements of ths structure*/
static rtwCAPI_DimensionMap rtDimensionMap[] = {
    /* dataOrientation, dimArrayIndex, numDims, vardimsIndex */
    {rtwCAPI_VECTOR, 0, 2, 0},

    {rtwCAPI_VECTOR, 2, 2, 0},

    {rtwCAPI_VECTOR, 4, 2, 0},

    {rtwCAPI_SCALAR, 6, 2, 0},

    {rtwCAPI_VECTOR, 8, 2, 0},

    {rtwCAPI_VECTOR, 10, 2, 0},

    {rtwCAPI_VECTOR, 12, 2, 0},

    {rtwCAPI_VECTOR, 14, 2, 0},

    {rtwCAPI_VECTOR, 16, 2, 0},

    {rtwCAPI_VECTOR, 18, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 20, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 22, 2, 0},

    {rtwCAPI_VECTOR, 24, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 26, 2, 0},

    {rtwCAPI_VECTOR, 28, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 30, 2, 0},

    {rtwCAPI_VECTOR, 32, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 34, 2, 0},

    {rtwCAPI_VECTOR, 36, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 38, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 40, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 42, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 6, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 44, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 46, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 48, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 50, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 52, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 10, 2, 0},

    {rtwCAPI_MATRIX_COL_MAJOR, 14, 2, 0}};

/* Dimension Array- use dimArrayIndex to access elements of this array */
static uint_T rtDimensionArray[] = {
    23, /* 0 */
    1, /* 1 */
    9, /* 2 */
    1, /* 3 */
    2, /* 4 */
    1, /* 5 */
    1, /* 6 */
    1, /* 7 */
    6, /* 8 */
    1, /* 9 */
    3, /* 10 */
    1, /* 11 */
    16, /* 12 */
    1, /* 13 */
    1, /* 14 */
    23, /* 15 */
    1, /* 16 */
    3, /* 17 */
    133, /* 18 */
    1, /* 19 */
    6, /* 20 */
    6, /* 21 */
    19, /* 22 */
    6, /* 23 */
    19, /* 24 */
    1, /* 25 */
    12, /* 26 */
    12, /* 27 */
    12, /* 28 */
    1, /* 29 */
    38, /* 30 */
    12, /* 31 */
    38, /* 32 */
    1, /* 33 */
    4, /* 34 */
    4, /* 35 */
    29, /* 36 */
    1, /* 37 */
    23, /* 38 */
    23, /* 39 */
    29, /* 40 */
    29, /* 41 */
    6, /* 42 */
    29, /* 43 */
    13, /* 44 */
    1, /* 45 */
    13, /* 46 */
    3, /* 47 */
    13, /* 48 */
    23, /* 49 */
    26, /* 50 */
    24, /* 51 */
    1, /* 52 */
    2 /* 53 */
};

/* C-API stores floating point values in an array. The elements of this  *
 * are unique. This ensures that values which are shared across the model*
 * are stored in the most efficient way. These values are referenced by  *
 *           - rtwCAPI_FixPtMap.fracSlopePtr,                            *
 *           - rtwCAPI_FixPtMap.biasPtr,                                 *
 *           - rtwCAPI_SampleTimeMap.samplePeriodPtr,                    *
 *           - rtwCAPI_SampleTimeMap.sampleOffsetPtr                     */
static const real_T rtcapiStoredFloats[] = {0.0, 0.01};

/* Fixed Point Map */
static rtwCAPI_FixPtMap rtFixPtMap[] = {
    /* fracSlopePtr, biasPtr, scaleType, wordLength, exponent, isSigned */
    {(NULL), (NULL), rtwCAPI_FIX_RESERVED, 0, 0, 0},
};

/* Sample Time Map - use sTimeIndex to access elements of ths structure */
static rtwCAPI_SampleTimeMap rtSampleTimeMap[] = {
    /* samplePeriodPtr, sampleOffsetPtr, tid, samplingMode */
    {(const void*) &rtcapiStoredFloats[0], (const void*) &rtcapiStoredFloats[0], 0, 0},

    {(const void*) &rtcapiStoredFloats[1], (const void*) &rtcapiStoredFloats[0], 1, 0}};

static rtwCAPI_ModelMappingStaticInfo mmiStatic = {
    /* Signals:{signals, numSignals,
     *           rootInputs, numRootInputs,
     *           rootOutputs, numRootOutputs},
     * Params: {blockParameters, numBlockParameters,
     *          modelParameters, numModelParameters},
     * States: {states, numStates},
     * Maps:   {dataTypeMap, dimensionMap, fixPtMap,
     *          elementMap, sampleTimeMap, dimensionArray},
     * TargetType: targetType
     */
    {rtBlockSignals, 157, (NULL), 0, (NULL), 0},

    {rtBlockParameters, 44, rtModelParameters, 8},

    {(NULL), 0},

    {rtDataTypeMap, rtDimensionMap, rtFixPtMap, rtElementMap, rtSampleTimeMap, rtDimensionArray},
    "float",

    {403027145U, 3969306854U, 2132273289U, 3949317986U},
    (NULL),
    0,
    0};

/* Function to get C API Model Mapping Static Info */
const rtwCAPI_ModelMappingStaticInfo* torqueBalancingYoga_GetCAPIStaticMap(void)
{
    return &mmiStatic;
}

/* Cache pointers into DataMapInfo substructure of RTModel */
#ifndef HOST_CAPI_BUILD

void torqueBalancingYoga_InitializeDataMapInfo(
    RT_MODEL_torqueBalancingYoga_T* const torqueBalancingYoga_M,
    B_torqueBalancingYoga_T* torqueBalancingYoga_B,
    P_torqueBalancingYoga_T* torqueBalancingYoga_P)
{
    /* Set C-API version */
    rtwCAPI_SetVersion(torqueBalancingYoga_M->DataMapInfo.mmi, 1);

    /* Cache static C-API data into the Real-time Model Data structure */
    rtwCAPI_SetStaticMap(torqueBalancingYoga_M->DataMapInfo.mmi, &mmiStatic);

    /* Cache static C-API logging data into the Real-time Model Data structure */
    rtwCAPI_SetLoggingStaticMap(torqueBalancingYoga_M->DataMapInfo.mmi, (NULL));

    /* Cache C-API Data Addresses into the Real-Time Model Data structure */
    torqueBalancingYoga_InitializeDataAddr(torqueBalancingYoga_M->DataMapInfo.dataAddress,
                                           torqueBalancingYoga_B,
                                           torqueBalancingYoga_P);
    rtwCAPI_SetDataAddressMap(torqueBalancingYoga_M->DataMapInfo.mmi,
                              torqueBalancingYoga_M->DataMapInfo.dataAddress);

    /* Cache C-API Data Run-Time Dimension Buffer Addresses into the Real-Time Model Data structure
     */
    torqueBalancingYoga_InitializeVarDimsAddr(torqueBalancingYoga_M->DataMapInfo.vardimsAddress);
    rtwCAPI_SetVarDimsAddressMap(torqueBalancingYoga_M->DataMapInfo.mmi,
                                 torqueBalancingYoga_M->DataMapInfo.vardimsAddress);

    /* Set Instance specific path */
    rtwCAPI_SetPath(torqueBalancingYoga_M->DataMapInfo.mmi, (NULL));
    rtwCAPI_SetFullPath(torqueBalancingYoga_M->DataMapInfo.mmi, (NULL));

    /* Cache C-API logging function pointers into the Real-Time Model Data structure */
    torqueBalancingYoga_InitializeLoggingFunctions(torqueBalancingYoga_M->DataMapInfo.loggingPtrs);
    rtwCAPI_SetLoggingPtrs(torqueBalancingYoga_M->DataMapInfo.mmi,
                           torqueBalancingYoga_M->DataMapInfo.loggingPtrs);

    /* Cache the instance C-API logging pointer */
    rtwCAPI_SetInstanceLoggingInfo(torqueBalancingYoga_M->DataMapInfo.mmi, (NULL));

    /* Set reference to submodels */
    rtwCAPI_SetChildMMIArray(torqueBalancingYoga_M->DataMapInfo.mmi, (NULL));
    rtwCAPI_SetChildMMIArrayLen(torqueBalancingYoga_M->DataMapInfo.mmi, 0);
}

#else /* HOST_CAPI_BUILD */
#ifdef __cplusplus

extern "C" {

#endif

void torqueBalancingYoga_host_InitializeDataMapInfo(torqueBalancingYoga_host_DataMapInfo_T* dataMap,
                                                    const char* path)
{
    /* Set C-API version */
    rtwCAPI_SetVersion(dataMap->mmi, 1);

    /* Cache static C-API data into the Real-time Model Data structure */
    rtwCAPI_SetStaticMap(dataMap->mmi, &mmiStatic);

    /* host data address map is NULL */
    rtwCAPI_SetDataAddressMap(dataMap->mmi, NULL);

    /* host vardims address map is NULL */
    rtwCAPI_SetVarDimsAddressMap(dataMap->mmi, NULL);

    /* Set Instance specific path */
    rtwCAPI_SetPath(dataMap->mmi, path);
    rtwCAPI_SetFullPath(dataMap->mmi, NULL);

    /* Set reference to submodels */
    rtwCAPI_SetChildMMIArray(dataMap->mmi, (NULL));
    rtwCAPI_SetChildMMIArrayLen(dataMap->mmi, 0);
}

#ifdef __cplusplus
}
#endif
#endif /* HOST_CAPI_BUILD */

/* EOF: torqueBalancingYoga_capi.cpp */
