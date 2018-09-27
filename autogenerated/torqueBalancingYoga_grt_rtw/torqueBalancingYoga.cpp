/*
 * torqueBalancingYoga.cpp
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueBalancingYoga".
 *
 * Model version              : 1.3292
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C++ source code generated on : Thu Sep 27 15:05:09 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "torqueBalancingYoga.h"
#include "torqueBalancingYoga_private.h"

/* Forward declaration for local functions */
static void torqueBalancingYoga_rotz(real_T alpha, real_T R[9]);
static void torqueBalancingYoga_roty(real_T alpha, real_T R[9]);
static void torqueBalancingYoga_xgetrf(real_T A[16], int32_T ipiv[4], int32_T* info);

/* Function for MATLAB Function: '<S41>/fromImuToHomogeousTransformFCN' */
static void torqueBalancingYoga_rotz(real_T alpha, real_T R[9])
{
    real_T R_tmp;
    real_T R_tmp_0;
    memset(&R[0], 0, 9U * sizeof(real_T));
    R[8] = 1.0;
    R_tmp_0 = std::cos(alpha);
    R[0] = R_tmp_0;
    R_tmp = std::sin(alpha);
    R[3] = -R_tmp;
    R[1] = R_tmp;
    R[4] = R_tmp_0;
}

/* Function for MATLAB Function: '<S41>/fromImuToHomogeousTransformFCN' */
static void torqueBalancingYoga_roty(real_T alpha, real_T R[9])
{
    real_T R_tmp;
    real_T R_tmp_0;
    memset(&R[0], 0, 9U * sizeof(real_T));
    R[4] = 1.0;
    R_tmp_0 = std::cos(alpha);
    R[0] = R_tmp_0;
    R_tmp = std::sin(alpha);
    R[6] = R_tmp;
    R[2] = -R_tmp;
    R[8] = R_tmp_0;
}

/* Function for MATLAB Function: '<S41>/fromImuToHomogeousTransformFCN' */
static void torqueBalancingYoga_xgetrf(real_T A[16], int32_T ipiv[4], int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    ipiv[0] = 1;
    ipiv[1] = 2;
    ipiv[2] = 3;
    ipiv[3] = 4;
    *info = 0;
    for (j = 0; j < 3; j++) {
        c = j * 5;
        iy = 0;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 4 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                iy = k - 1;
                smax = s;
            }
        }

        if (A[c + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = iy + 1;
                smax = A[j];
                A[j] = A[iy];
                A[iy] = smax;
                ix = j + 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
                ix += 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
                ix += 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
            }

            iy = (c - j) + 4;
            for (ix = c + 1; ix < iy; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        iy = c;
        ix = c + 4;
        for (k = 1; k <= 3 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                c_ix = c + 1;
                d = (iy - j) + 8;
                for (ijA = 5 + iy; ijA < d; ijA++) {
                    A[ijA] += A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 4;
            iy += 4;
        }
    }

    if ((*info == 0) && (!(A[15] != 0.0))) {
        *info = 4;
    }
}

/*
 * Output and update for atomic system:
 *    '<S41>/fromImuToHomogeousTransformFCN'
 *    '<S78>/fromImuToHomogeousTransformFCN'
 *    '<S80>/fromImuToHomogeousTransformFCN'
 *    '<S118>/fromImuToHomogeousTransformFCN'
 *    '<S120>/fromImuToHomogeousTransformFCN'
 */
void fromImuToHomogeousTransformFCN(const real_T rtu_imu_H_link[16],
                                    const real_T rtu_imu_H_link_0[16],
                                    const real_T rtu_link_H_base[16],
                                    const real_T rtu_inertial_0[12],
                                    const real_T rtu_inertial[12],
                                    const real_T rtu_neck_pos[3],
                                    B_fromImuToHomogeousTransform_T* localB,
                                    P_torqueBalancingYoga_T* torqueBalancingYoga_P)
{
    real_T inertial_0[12];
    real_T inertial[12];
    real_T wImu_R_link[9];
    real_T wImu_R_link_0[9];
    real_T rollPitchYaw_link_0[3];
    real_T wImu_H_wImuAssumingNeckToZero[16];
    real_T R[9];
    real_T b_R[9];
    real_T c_R[9];
    real_T d_R[9];
    real_T A[16];
    int32_T ipiv[4];
    int32_T info;
    int32_T jAcol;
    int32_T kBcol;
    int32_T kAcol;
    static const real_T B[16] = {-1.0,
                                 1.2246467991473532E-16,
                                 -1.2246467991473532E-16,
                                 0.0,
                                 -1.2246467991473532E-16,
                                 -6.1232339957367648E-17,
                                 1.0,
                                 0.0,
                                 1.2246467991473532E-16,
                                 1.0,
                                 6.123233995736766E-17,
                                 0.0,
                                 -0.018499999999999985,
                                 0.12029999999999999,
                                 0.0066000000000000043,
                                 1.0};

    static const real_T b[16] = {1.0,
                                 0.0,
                                 0.0,
                                 0.0,
                                 -0.0,
                                 6.123233995736766E-17,
                                 1.0,
                                 0.0,
                                 0.0,
                                 -1.0,
                                 6.123233995736766E-17,
                                 0.0,
                                 0.0,
                                 0.0,
                                 0.0066,
                                 1.0};

    static const int8_T c[4] = {0, 0, 0, 1};

    static const real_T d[4] = {0.0, -1.0, 6.123233995736766E-17, 0.1108};

    int32_T i;
    real_T tmp[9];
    real_T tmp_0[9];
    real_T tmp_1[9];
    real_T b_R_0[9];
    real_T b_R_1[9];
    real_T tmp_2[16];
    real_T tmp_3[16];
    real_T R_0[9];
    real_T R_1[16];
    real_T R_2[16];
    real_T rollPitchYaw_link_idx_0;
    real_T rollPitchYawFiltered_link_idx_2;
    real_T rollPitchYawFiltered_link_idx_1;
    real_T R_tmp;

    /* MATLAB Function 'Robot State and References/Internal Coordinator/Base to fixed_link/LFoot to
     * base link transform /fromImuToHomogeousTransformFCN': '<S46>:1' */
    /* '<S46>:1:3' */
    for (i = 0; i < 12; i++) {
        inertial[i] = rtu_inertial[i] * 3.1415926535897931 / 180.0;
        inertial_0[i] = rtu_inertial_0[i] * 3.1415926535897931 / 180.0;
    }

    memset(&R[0], 0, 9U * sizeof(real_T));
    memset(&b_R[0], 0, 9U * sizeof(real_T));
    memset(&c_R[0], 0, 9U * sizeof(real_T));
    memset(&d_R[0], 0, 9U * sizeof(real_T));
    R[0] = 1.0;
    R_tmp = std::cos(inertial[0]);
    R[4] = R_tmp;
    rollPitchYawFiltered_link_idx_2 = std::sin(inertial[0]);
    R[7] = -rollPitchYawFiltered_link_idx_2;
    R[5] = rollPitchYawFiltered_link_idx_2;
    R[8] = R_tmp;
    b_R[8] = 1.0;
    R_tmp = std::cos(inertial_0[2]);
    b_R[0] = R_tmp;
    rollPitchYawFiltered_link_idx_1 = std::sin(inertial_0[2]);
    b_R[3] = -rollPitchYawFiltered_link_idx_1;
    b_R[1] = rollPitchYawFiltered_link_idx_1;
    b_R[4] = R_tmp;
    c_R[4] = 1.0;
    rollPitchYawFiltered_link_idx_1 = std::cos(inertial_0[1]);
    c_R[0] = rollPitchYawFiltered_link_idx_1;
    rollPitchYaw_link_idx_0 = std::sin(inertial_0[1]);
    c_R[6] = rollPitchYaw_link_idx_0;
    c_R[2] = -rollPitchYaw_link_idx_0;
    c_R[8] = rollPitchYawFiltered_link_idx_1;
    d_R[0] = 1.0;
    rollPitchYaw_link_idx_0 = std::cos(inertial_0[0]);
    d_R[4] = rollPitchYaw_link_idx_0;
    rollPitchYawFiltered_link_idx_1 = std::sin(inertial_0[0]);
    d_R[7] = -rollPitchYawFiltered_link_idx_1;
    d_R[5] = rollPitchYawFiltered_link_idx_1;
    d_R[8] = rollPitchYaw_link_idx_0;
    torqueBalancingYoga_rotz(inertial[2], R_0);
    torqueBalancingYoga_roty(inertial[1], tmp);
    for (i = 0; i < 3; i++) {
        for (info = 0; info < 3; info++) {
            jAcol = i + 3 * info;
            tmp_0[jAcol] = 0.0;
            tmp_0[jAcol] = tmp_0[3 * info + i] + tmp[3 * info] * R_0[i];
            tmp_0[jAcol] = tmp[3 * info + 1] * R_0[i + 3] + tmp_0[3 * info + i];
            tmp_0[jAcol] = tmp[3 * info + 2] * R_0[i + 6] + tmp_0[3 * info + i];
        }

        for (info = 0; info < 3; info++) {
            jAcol = i + 3 * info;
            tmp_1[jAcol] = 0.0;
            tmp_1[jAcol] = tmp_1[3 * info + i] + R[3 * info] * tmp_0[i];
            tmp_1[jAcol] = R[3 * info + 1] * tmp_0[i + 3] + tmp_1[3 * info + i];
            tmp_1[jAcol] = R[3 * info + 2] * tmp_0[i + 6] + tmp_1[3 * info + i];
        }

        for (info = 0; info < 3; info++) {
            jAcol = i + 3 * info;
            wImu_R_link[jAcol] = 0.0;
            b_R_0[jAcol] = 0.0;
            kAcol = 3 * info + i;
            wImu_R_link[jAcol] = wImu_R_link[kAcol] + rtu_imu_H_link[info << 2] * tmp_1[i];
            b_R_0[jAcol] = b_R_0[kAcol] + c_R[3 * info] * b_R[i];
            wImu_R_link[jAcol] =
                rtu_imu_H_link[(info << 2) + 1] * tmp_1[i + 3] + wImu_R_link[3 * info + i];
            b_R_0[jAcol] = c_R[3 * info + 1] * b_R[i + 3] + b_R_0[3 * info + i];
            wImu_R_link[jAcol] =
                rtu_imu_H_link[(info << 2) + 2] * tmp_1[i + 6] + wImu_R_link[3 * info + i];
            b_R_0[jAcol] = c_R[3 * info + 2] * b_R[i + 6] + b_R_0[3 * info + i];
        }

        for (info = 0; info < 3; info++) {
            jAcol = i + 3 * info;
            b_R_1[jAcol] = 0.0;
            b_R_1[jAcol] = b_R_1[3 * info + i] + d_R[3 * info] * b_R_0[i];
            b_R_1[jAcol] = d_R[3 * info + 1] * b_R_0[i + 3] + b_R_1[3 * info + i];
            b_R_1[jAcol] = d_R[3 * info + 2] * b_R_0[i + 6] + b_R_1[3 * info + i];
        }

        for (info = 0; info < 3; info++) {
            jAcol = i + 3 * info;
            wImu_R_link_0[jAcol] = 0.0;
            wImu_R_link_0[jAcol] =
                wImu_R_link_0[3 * info + i] + rtu_imu_H_link_0[info << 2] * b_R_1[i];
            wImu_R_link_0[jAcol] =
                rtu_imu_H_link_0[(info << 2) + 1] * b_R_1[i + 3] + wImu_R_link_0[3 * info + i];
            wImu_R_link_0[jAcol] =
                rtu_imu_H_link_0[(info << 2) + 2] * b_R_1[i + 6] + wImu_R_link_0[3 * info + i];
        }

        rollPitchYaw_link_0[i] = 0.0;
    }

    if (wImu_R_link_0[2] < 1.0) {
        if (wImu_R_link_0[2] > -1.0) {
            rollPitchYaw_link_0[1] = std::asin(-wImu_R_link_0[2]);
            rollPitchYaw_link_0[2] = atan2(wImu_R_link_0[1], wImu_R_link_0[0]);
        }
        else {
            rollPitchYaw_link_0[2] = -atan2(-wImu_R_link_0[7], wImu_R_link_0[4]);
        }
    }
    else {
        rollPitchYaw_link_0[2] = atan2(-wImu_R_link_0[7], wImu_R_link_0[4]);
    }

    rollPitchYaw_link_idx_0 = 0.0;
    rollPitchYawFiltered_link_idx_1 = 0.0;
    if (wImu_R_link[2] < 1.0) {
        if (wImu_R_link[2] > -1.0) {
            rollPitchYawFiltered_link_idx_1 = std::asin(-wImu_R_link[2]);
            rollPitchYawFiltered_link_idx_2 = atan2(wImu_R_link[1], wImu_R_link[0]);
            rollPitchYaw_link_idx_0 = atan2(wImu_R_link[5], wImu_R_link[8]);
        }
        else {
            rollPitchYawFiltered_link_idx_2 = -atan2(-wImu_R_link[7], wImu_R_link[4]);
        }
    }
    else {
        rollPitchYawFiltered_link_idx_2 = atan2(-wImu_R_link[7], wImu_R_link[4]);
    }

    if (torqueBalancingYoga_P->Config.FILTER_IMU_YAW) {
        rollPitchYawFiltered_link_idx_2 = rollPitchYaw_link_0[2];
    }

    if (torqueBalancingYoga_P->Config.FILTER_IMU_PITCH) {
        rollPitchYawFiltered_link_idx_1 = rollPitchYaw_link_0[1];
    }

    memset(&R[0], 0, 9U * sizeof(real_T));
    memset(&b_R[0], 0, 9U * sizeof(real_T));
    memset(&c_R[0], 0, 9U * sizeof(real_T));
    R[8] = 1.0;
    R_tmp = std::cos(rollPitchYawFiltered_link_idx_2);
    R[0] = R_tmp;
    rollPitchYawFiltered_link_idx_2 = std::sin(rollPitchYawFiltered_link_idx_2);
    R[3] = -rollPitchYawFiltered_link_idx_2;
    R[1] = rollPitchYawFiltered_link_idx_2;
    R[4] = R_tmp;
    b_R[4] = 1.0;
    R_tmp = std::cos(rollPitchYawFiltered_link_idx_1);
    b_R[0] = R_tmp;
    rollPitchYawFiltered_link_idx_1 = std::sin(rollPitchYawFiltered_link_idx_1);
    b_R[6] = rollPitchYawFiltered_link_idx_1;
    b_R[2] = -rollPitchYawFiltered_link_idx_1;
    b_R[8] = R_tmp;
    c_R[0] = 1.0;
    rollPitchYawFiltered_link_idx_1 = std::cos(rollPitchYaw_link_idx_0);
    c_R[4] = rollPitchYawFiltered_link_idx_1;
    rollPitchYaw_link_idx_0 = std::sin(rollPitchYaw_link_idx_0);
    c_R[7] = -rollPitchYaw_link_idx_0;
    c_R[5] = rollPitchYaw_link_idx_0;
    c_R[8] = rollPitchYawFiltered_link_idx_1;
    memcpy(&A[0], &B[0], sizeof(real_T) << 4U);
    torqueBalancingYoga_xgetrf(A, ipiv, &info);
    rollPitchYaw_link_idx_0 = std::cos(rtu_neck_pos[0] + 1.5707963267948966);
    R_1[0] = rollPitchYaw_link_idx_0;
    rollPitchYawFiltered_link_idx_1 = std::sin(rtu_neck_pos[0] + 1.5707963267948966);
    R_1[4] = -rollPitchYawFiltered_link_idx_1 * 6.123233995736766E-17;
    R_1[8] = rollPitchYawFiltered_link_idx_1;
    R_1[12] = rollPitchYaw_link_idx_0 * 0.0095;
    R_1[1] = rollPitchYawFiltered_link_idx_1;
    R_1[5] = rollPitchYaw_link_idx_0 * 6.123233995736766E-17;
    R_1[9] = -rollPitchYaw_link_idx_0;
    R_1[13] = rollPitchYawFiltered_link_idx_1 * 0.0095;
    rollPitchYaw_link_idx_0 = std::cos(rtu_neck_pos[1] - 1.5707963267948966);
    R_2[0] = rollPitchYaw_link_idx_0;
    rollPitchYawFiltered_link_idx_1 = std::sin(rtu_neck_pos[1] - 1.5707963267948966);
    R_2[4] = -rollPitchYawFiltered_link_idx_1 * 6.123233995736766E-17;
    R_2[8] = -rollPitchYawFiltered_link_idx_1;
    R_2[12] = 0.0;
    R_2[1] = rollPitchYawFiltered_link_idx_1;
    R_2[5] = rollPitchYaw_link_idx_0 * 6.123233995736766E-17;
    R_2[9] = -(-rollPitchYaw_link_idx_0);
    R_2[13] = 0.0;
    R_1[2] = 0.0;
    R_1[3] = 0.0;
    R_2[2] = 0.0;
    R_2[3] = 0.0;
    R_1[6] = 1.0;
    R_1[7] = 0.0;
    R_2[6] = -1.0;
    R_2[7] = 0.0;
    R_1[10] = 6.123233995736766E-17;
    R_1[11] = 0.0;
    R_2[10] = 6.123233995736766E-17;
    R_2[11] = 0.0;
    R_1[14] = 0.0;
    R_1[15] = 1.0;
    R_2[14] = 0.0;
    R_2[15] = 1.0;
    rollPitchYaw_link_idx_0 = std::cos(rtu_neck_pos[2] + 1.5707963267948966);
    tmp_3[0] = rollPitchYaw_link_idx_0;
    rollPitchYawFiltered_link_idx_1 = std::sin(rtu_neck_pos[2] + 1.5707963267948966);
    tmp_3[4] = -rollPitchYawFiltered_link_idx_1 * 6.123233995736766E-17;
    tmp_3[8] = -rollPitchYawFiltered_link_idx_1;
    tmp_3[12] = rollPitchYaw_link_idx_0 * 0.0185;
    tmp_3[1] = rollPitchYawFiltered_link_idx_1;
    tmp_3[5] = rollPitchYaw_link_idx_0 * 6.123233995736766E-17;
    tmp_3[9] = -(-rollPitchYaw_link_idx_0);
    tmp_3[13] = rollPitchYawFiltered_link_idx_1 * 0.0185;
    for (info = 0; info < 4; info++) {
        for (jAcol = 0; jAcol < 4; jAcol++) {
            tmp_2[info + (jAcol << 2)] = 0.0;
            tmp_2[info + (jAcol << 2)] += R_2[jAcol << 2] * R_1[info];
            tmp_2[info + (jAcol << 2)] += R_2[(jAcol << 2) + 1] * R_1[info + 4];
            tmp_2[info + (jAcol << 2)] += R_2[(jAcol << 2) + 2] * R_1[info + 8];
            tmp_2[info + (jAcol << 2)] += R_2[(jAcol << 2) + 3] * R_1[info + 12];
        }

        tmp_3[2 + (info << 2)] = d[info];
        tmp_3[3 + (info << 2)] = c[info];
    }

    for (info = 0; info < 4; info++) {
        for (jAcol = 0; jAcol < 4; jAcol++) {
            R_1[info + (jAcol << 2)] = 0.0;
            R_1[info + (jAcol << 2)] += tmp_3[jAcol << 2] * tmp_2[info];
            R_1[info + (jAcol << 2)] += tmp_3[(jAcol << 2) + 1] * tmp_2[info + 4];
            R_1[info + (jAcol << 2)] += tmp_3[(jAcol << 2) + 2] * tmp_2[info + 8];
            R_1[info + (jAcol << 2)] += tmp_3[(jAcol << 2) + 3] * tmp_2[info + 12];
        }

        for (jAcol = 0; jAcol < 4; jAcol++) {
            wImu_H_wImuAssumingNeckToZero[info + (jAcol << 2)] = 0.0;
            wImu_H_wImuAssumingNeckToZero[info + (jAcol << 2)] += b[jAcol << 2] * R_1[info];
            wImu_H_wImuAssumingNeckToZero[info + (jAcol << 2)] +=
                b[(jAcol << 2) + 1] * R_1[info + 4];
            wImu_H_wImuAssumingNeckToZero[info + (jAcol << 2)] +=
                b[(jAcol << 2) + 2] * R_1[info + 8];
            wImu_H_wImuAssumingNeckToZero[info + (jAcol << 2)] +=
                b[(jAcol << 2) + 3] * R_1[info + 12];
        }
    }

    for (info = 0; info < 4; info++) {
        i = info << 2;
        jAcol = info << 2;
        for (kAcol = 1; kAcol <= info; kAcol++) {
            kBcol = (kAcol - 1) << 2;
            if (A[(kAcol + jAcol) - 1] != 0.0) {
                wImu_H_wImuAssumingNeckToZero[i] -=
                    A[(kAcol + jAcol) - 1] * wImu_H_wImuAssumingNeckToZero[kBcol];
                wImu_H_wImuAssumingNeckToZero[1 + i] -=
                    A[(kAcol + jAcol) - 1] * wImu_H_wImuAssumingNeckToZero[1 + kBcol];
                wImu_H_wImuAssumingNeckToZero[2 + i] -=
                    A[(kAcol + jAcol) - 1] * wImu_H_wImuAssumingNeckToZero[2 + kBcol];
                wImu_H_wImuAssumingNeckToZero[3 + i] -=
                    A[(kAcol + jAcol) - 1] * wImu_H_wImuAssumingNeckToZero[3 + kBcol];
            }
        }

        rollPitchYaw_link_idx_0 = 1.0 / A[info + jAcol];
        wImu_H_wImuAssumingNeckToZero[i] *= rollPitchYaw_link_idx_0;
        wImu_H_wImuAssumingNeckToZero[1 + i] *= rollPitchYaw_link_idx_0;
        wImu_H_wImuAssumingNeckToZero[2 + i] *= rollPitchYaw_link_idx_0;
        wImu_H_wImuAssumingNeckToZero[3 + i] *= rollPitchYaw_link_idx_0;
    }

    for (info = 3; info >= 0; info--) {
        i = info << 2;
        jAcol = (info << 2) - 1;
        for (kAcol = info + 2; kAcol < 5; kAcol++) {
            kBcol = (kAcol - 1) << 2;
            if (A[kAcol + jAcol] != 0.0) {
                wImu_H_wImuAssumingNeckToZero[i] -=
                    A[kAcol + jAcol] * wImu_H_wImuAssumingNeckToZero[kBcol];
                wImu_H_wImuAssumingNeckToZero[1 + i] -=
                    A[kAcol + jAcol] * wImu_H_wImuAssumingNeckToZero[1 + kBcol];
                wImu_H_wImuAssumingNeckToZero[2 + i] -=
                    A[kAcol + jAcol] * wImu_H_wImuAssumingNeckToZero[2 + kBcol];
                wImu_H_wImuAssumingNeckToZero[3 + i] -=
                    A[kAcol + jAcol] * wImu_H_wImuAssumingNeckToZero[3 + kBcol];
            }
        }
    }

    for (info = 2; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            i = ipiv[info] - 1;
            rollPitchYaw_link_idx_0 = wImu_H_wImuAssumingNeckToZero[info << 2];
            wImu_H_wImuAssumingNeckToZero[info << 2] = wImu_H_wImuAssumingNeckToZero[i << 2];
            wImu_H_wImuAssumingNeckToZero[i << 2] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = wImu_H_wImuAssumingNeckToZero[(info << 2) + 1];
            wImu_H_wImuAssumingNeckToZero[1 + (info << 2)] =
                wImu_H_wImuAssumingNeckToZero[(i << 2) + 1];
            wImu_H_wImuAssumingNeckToZero[1 + (i << 2)] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = wImu_H_wImuAssumingNeckToZero[(info << 2) + 2];
            wImu_H_wImuAssumingNeckToZero[2 + (info << 2)] =
                wImu_H_wImuAssumingNeckToZero[(i << 2) + 2];
            wImu_H_wImuAssumingNeckToZero[2 + (i << 2)] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = wImu_H_wImuAssumingNeckToZero[(info << 2) + 3];
            wImu_H_wImuAssumingNeckToZero[3 + (info << 2)] =
                wImu_H_wImuAssumingNeckToZero[(i << 2) + 3];
            wImu_H_wImuAssumingNeckToZero[3 + (i << 2)] = rollPitchYaw_link_idx_0;
        }
    }

    for (info = 0; info < 3; info++) {
        A[info << 2] = wImu_R_link_0[3 * info];
        A[1 + (info << 2)] = wImu_R_link_0[3 * info + 1];
        A[2 + (info << 2)] = wImu_R_link_0[3 * info + 2];
        A[12 + info] = 0.0;
    }

    A[3] = 0.0;
    A[7] = 0.0;
    A[11] = 0.0;
    A[15] = 1.0;
    torqueBalancingYoga_xgetrf(A, ipiv, &info);
    for (info = 0; info < 3; info++) {
        for (jAcol = 0; jAcol < 3; jAcol++) {
            i = info + 3 * jAcol;
            d_R[i] = 0.0;
            d_R[i] = d_R[3 * jAcol + info] + b_R[3 * jAcol] * R[info];
            d_R[i] = b_R[3 * jAcol + 1] * R[info + 3] + d_R[3 * jAcol + info];
            d_R[i] = b_R[3 * jAcol + 2] * R[info + 6] + d_R[3 * jAcol + info];
        }

        for (jAcol = 0; jAcol < 3; jAcol++) {
            i = info + 3 * jAcol;
            R_0[i] = 0.0;
            R_0[i] = R_0[3 * jAcol + info] + c_R[3 * jAcol] * d_R[info];
            R_0[i] = c_R[3 * jAcol + 1] * d_R[info + 3] + R_0[3 * jAcol + info];
            R_0[i] = c_R[3 * jAcol + 2] * d_R[info + 6] + R_0[3 * jAcol + info];
        }
    }

    for (info = 0; info < 3; info++) {
        R_1[info << 2] = R_0[3 * info];
        R_1[1 + (info << 2)] = R_0[3 * info + 1];
        R_1[2 + (info << 2)] = R_0[3 * info + 2];
        R_1[12 + info] = 0.0;
    }

    R_1[3] = 0.0;
    R_1[7] = 0.0;
    R_1[11] = 0.0;
    R_1[15] = 1.0;
    for (info = 0; info < 4; info++) {
        for (jAcol = 0; jAcol < 4; jAcol++) {
            R_2[jAcol + (info << 2)] = 0.0;
            R_2[jAcol + (info << 2)] += rtu_link_H_base[info << 2] * R_1[jAcol];
            R_2[jAcol + (info << 2)] += rtu_link_H_base[(info << 2) + 1] * R_1[jAcol + 4];
            R_2[jAcol + (info << 2)] += rtu_link_H_base[(info << 2) + 2] * R_1[jAcol + 8];
            R_2[jAcol + (info << 2)] += rtu_link_H_base[(info << 2) + 3] * R_1[jAcol + 12];
        }
    }

    for (info = 0; info < 4; info++) {
        for (jAcol = 0; jAcol < 4; jAcol++) {
            localB->w_H_b[jAcol + (info << 2)] = 0.0;
            localB->w_H_b[jAcol + (info << 2)] +=
                R_2[info << 2] * wImu_H_wImuAssumingNeckToZero[jAcol];
            localB->w_H_b[jAcol + (info << 2)] +=
                R_2[(info << 2) + 1] * wImu_H_wImuAssumingNeckToZero[jAcol + 4];
            localB->w_H_b[jAcol + (info << 2)] +=
                R_2[(info << 2) + 2] * wImu_H_wImuAssumingNeckToZero[jAcol + 8];
            localB->w_H_b[jAcol + (info << 2)] +=
                R_2[(info << 2) + 3] * wImu_H_wImuAssumingNeckToZero[jAcol + 12];
        }
    }

    for (info = 0; info < 3; info++) {
        if (info + 1 != ipiv[info]) {
            i = ipiv[info] - 1;
            rollPitchYaw_link_idx_0 = localB->w_H_b[info];
            localB->w_H_b[info] = localB->w_H_b[i];
            localB->w_H_b[i] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = localB->w_H_b[info + 4];
            localB->w_H_b[info + 4] = localB->w_H_b[i + 4];
            localB->w_H_b[i + 4] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = localB->w_H_b[info + 8];
            localB->w_H_b[info + 8] = localB->w_H_b[i + 8];
            localB->w_H_b[i + 8] = rollPitchYaw_link_idx_0;
            rollPitchYaw_link_idx_0 = localB->w_H_b[info + 12];
            localB->w_H_b[info + 12] = localB->w_H_b[i + 12];
            localB->w_H_b[i + 12] = rollPitchYaw_link_idx_0;
        }
    }

    for (info = 0; info < 4; info++) {
        i = info << 2;
        if (localB->w_H_b[i] != 0.0) {
            for (kBcol = 1; kBcol + 1 < 5; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[i] * A[kBcol];
            }
        }

        if (localB->w_H_b[1 + i] != 0.0) {
            for (kBcol = 2; kBcol + 1 < 5; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[1 + i] * A[kBcol + 4];
            }
        }

        if (localB->w_H_b[2 + i] != 0.0) {
            for (kBcol = 3; kBcol + 1 < 5; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[2 + i] * A[kBcol + 8];
            }
        }
    }

    for (info = 0; info < 4; info++) {
        i = info << 2;
        if (localB->w_H_b[3 + i] != 0.0) {
            localB->w_H_b[3 + i] /= A[15];
            for (kBcol = 0; kBcol < 3; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[3 + i] * A[kBcol + 12];
            }
        }

        if (localB->w_H_b[2 + i] != 0.0) {
            localB->w_H_b[2 + i] /= A[10];
            for (kBcol = 0; kBcol < 2; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[2 + i] * A[kBcol + 8];
            }
        }

        if (localB->w_H_b[1 + i] != 0.0) {
            localB->w_H_b[1 + i] /= A[5];
            for (kBcol = 0; kBcol < 1; kBcol++) {
                localB->w_H_b[kBcol + i] -= localB->w_H_b[1 + i] * A[kBcol + 4];
            }
        }

        if (localB->w_H_b[i] != 0.0) {
            localB->w_H_b[i] /= A[0];
        }
    }
}

/*
 * System initialize for atomic system:
 *    '<S47>/MATLAB Function'
 *    '<S85>/MATLAB Function'
 *    '<S95>/MATLAB Function'
 *    '<S125>/MATLAB Function'
 *    '<S135>/MATLAB Function'
 */
void torqueBalan_MATLABFunction_Init(DW_MATLABFunction_torqueBalan_T* localDW)
{
    localDW->state_not_empty = false;
}

/*
 * Output and update for atomic system:
 *    '<S47>/MATLAB Function'
 *    '<S85>/MATLAB Function'
 *    '<S95>/MATLAB Function'
 *    '<S125>/MATLAB Function'
 *    '<S135>/MATLAB Function'
 */
void torqueBalancingY_MATLABFunction(const real_T rtu_s[16],
                                     B_MATLABFunction_torqueBalanc_T* localB,
                                     DW_MATLABFunction_torqueBalan_T* localDW)
{
    /* MATLAB Function 'Utilities/holder /MATLAB Function': '<S51>:1' */
    if (!localDW->state_not_empty) {
        /* '<S51>:1:5' */
        /* '<S51>:1:6' */
        memcpy(&localDW->state[0], &rtu_s[0], sizeof(real_T) << 4U);
        localDW->state_not_empty = true;
    }

    /* '<S51>:1:9' */
    memcpy(&localB->s0[0], &localDW->state[0], sizeof(real_T) << 4U);
}

/*
 * System initialize for atomic system:
 *    '<S48>/MATLAB Function'
 *    '<S86>/MATLAB Function'
 *    '<S96>/MATLAB Function'
 *    '<S126>/MATLAB Function'
 *    '<S136>/MATLAB Function'
 */
void torqueBal_MATLABFunction_m_Init(DW_MATLABFunction_torqueBal_k_T* localDW)
{
    localDW->state_not_empty = false;
}

/*
 * Output and update for atomic system:
 *    '<S48>/MATLAB Function'
 *    '<S86>/MATLAB Function'
 *    '<S96>/MATLAB Function'
 *    '<S126>/MATLAB Function'
 *    '<S136>/MATLAB Function'
 */
void torqueBalancin_MATLABFunction_l(const real_T rtu_s[12],
                                     B_MATLABFunction_torqueBala_f_T* localB,
                                     DW_MATLABFunction_torqueBal_k_T* localDW)
{
    /* MATLAB Function 'Utilities/holder /MATLAB Function': '<S53>:1' */
    if (!localDW->state_not_empty) {
        /* '<S53>:1:5' */
        /* '<S53>:1:6' */
        memcpy(&localDW->state[0], &rtu_s[0], 12U * sizeof(real_T));
        localDW->state_not_empty = true;
    }

    /* '<S53>:1:9' */
    memcpy(&localB->s0[0], &localDW->state[0], 12U * sizeof(real_T));
}

/*
 * System initialize for atomic system:
 *    '<S38>/MATLAB Function'
 *    '<S69>/MATLAB Function'
 *    '<S162>/MATLAB Function'
 */
void torqueBa_MATLABFunction_mb_Init(DW_MATLABFunction_torqueBal_l_T* localDW)
{
    localDW->state_not_empty = false;
}

/*
 * Output and update for atomic system:
 *    '<S38>/MATLAB Function'
 *    '<S69>/MATLAB Function'
 *    '<S162>/MATLAB Function'
 */
void torqueBalancin_MATLABFunction_k(const real_T rtu_s[23],
                                     B_MATLABFunction_torqueBala_j_T* localB,
                                     DW_MATLABFunction_torqueBal_l_T* localDW)
{
    /* MATLAB Function 'Utilities/holder /MATLAB Function': '<S60>:1' */
    if (!localDW->state_not_empty) {
        /* '<S60>:1:5' */
        /* '<S60>:1:6' */
        memcpy(&localDW->state[0], &rtu_s[0], 23U * sizeof(real_T));
        localDW->state_not_empty = true;
    }

    /* '<S60>:1:9' */
    memcpy(&localB->s0[0], &localDW->state[0], 23U * sizeof(real_T));
}

/*
 * System initialize for atomic system:
 *    '<S39>/MATLAB Function'
 *    '<S68>/MATLAB Function'
 */
void torqueBal_MATLABFunction_l_Init(DW_MATLABFunction_torqueBal_e_T* localDW)
{
    localDW->state_not_empty = false;
}

/*
 * Output and update for atomic system:
 *    '<S39>/MATLAB Function'
 *    '<S68>/MATLAB Function'
 */
void torqueBalancin_MATLABFunction_a(const real_T rtu_s[3],
                                     B_MATLABFunction_torqueBala_p_T* localB,
                                     DW_MATLABFunction_torqueBal_e_T* localDW)
{
    /* MATLAB Function 'Utilities/holder /MATLAB Function': '<S62>:1' */
    if (!localDW->state_not_empty) {
        /* '<S62>:1:5' */
        /* '<S62>:1:6' */
        localDW->state[0] = rtu_s[0];
        localDW->state[1] = rtu_s[1];
        localDW->state[2] = rtu_s[2];
        localDW->state_not_empty = true;
    }

    /* '<S62>:1:9' */
    localB->s0[0] = localDW->state[0];
    localB->s0[1] = localDW->state[1];
    localB->s0[2] = localDW->state[2];
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xswap(real_T x[529],
                                                              int32_T ix0,
                                                              int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 23; k++) {
        temp = x[ix];
        x[ix] = x[iy];
        x[iy] = temp;
        ix += 23;
        iy += 23;
    }
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xgetrf_h(real_T A[529],
                                                                 int32_T ipiv[23],
                                                                 int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T jA;
    int32_T b_ix;
    int32_T d;
    int32_T ijA;
    for (j = 0; j < 23; j++) {
        ipiv[j] = 1 + j;
    }

    *info = 0;
    for (j = 0; j < 22; j++) {
        c = j * 24;
        jA = 1;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 23 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                jA = k;
                smax = s;
            }
        }

        if (A[(c + jA) - 1] != 0.0) {
            if (jA - 1 != 0) {
                ipiv[j] = j + jA;
                torqueBalancingYoga_xswap(A, j + 1, j + jA);
            }

            jA = (c - j) + 23;
            for (ix = c + 1; ix < jA; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        jA = c;
        ix = c + 23;
        for (k = 1; k <= 22 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                b_ix = c + 1;
                d = (jA - j) + 46;
                for (ijA = 24 + jA; ijA < d; ijA++) {
                    A[ijA] += A[b_ix] * -smax;
                    b_ix++;
                }
            }

            ix += 23;
            jA += 23;
        }
    }

    if ((*info == 0) && (!(A[528] != 0.0))) {
        *info = 23;
    }
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm(const real_T A[529], real_T B[529])
{
    real_T temp;
    int32_T jBcol;
    int32_T jAcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 0; j < 23; j++) {
        jBcol = 23 * j;
        jAcol = 23 * j;
        for (k = 1; k <= j; k++) {
            kBcol = (k - 1) * 23;
            if (A[(k + jAcol) - 1] != 0.0) {
                for (i = 0; i < 23; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[(k + jAcol) - 1] * B[i + kBcol];
                }
            }
        }

        temp = 1.0 / A[j + jAcol];
        for (jAcol = 0; jAcol < 23; jAcol++) {
            tmp = jAcol + jBcol;
            B[tmp] *= temp;
        }
    }
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm_h(const real_T A[529], real_T B[529])
{
    int32_T jAcol;
    int32_T jBcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 22; j >= 0; j--) {
        jBcol = 23 * j;
        jAcol = 23 * j - 1;
        for (k = j + 2; k < 24; k++) {
            kBcol = (k - 1) * 23;
            if (A[k + jAcol] != 0.0) {
                for (i = 0; i < 23; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[k + jAcol] * B[i + kBcol];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mrdivide(const real_T A[529],
                                                                 const real_T B[529],
                                                                 real_T y[529])
{
    int32_T jp;
    int32_T ipiv[23];
    int32_T info;
    real_T temp;
    int32_T xi;
    memcpy(&torqueBalancingYoga_B.b_A_k[0], &B[0], 529U * sizeof(real_T));
    torqueBalancingYoga_xgetrf_h(torqueBalancingYoga_B.b_A_k, ipiv, &info);
    memcpy(&y[0], &A[0], 529U * sizeof(real_T));
    torqueBalancingYoga_xtrsm(torqueBalancingYoga_B.b_A_k, y);
    torqueBalancingYoga_xtrsm_h(torqueBalancingYoga_B.b_A_k, y);
    for (info = 21; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            jp = ipiv[info] - 1;
            for (xi = 0; xi < 23; xi++) {
                temp = y[23 * info + xi];
                y[xi + 23 * info] = y[23 * jp + xi];
                y[xi + 23 * jp] = temp;
            }
        }
    }
}

void rt_mldivided4x4(const real_T u0[16], const real_T u1[16], real_T y[16])
{
    real_T temp;
    real_T A[16];
    int8_T ipiv[4];
    int32_T j;
    int32_T ix;
    real_T s;
    int32_T iy;
    int32_T j_0;
    int32_T ijA;
    int32_T jBcol;
    int32_T kAcol;
    int32_T k;
    memcpy(&A[0], &u0[0], sizeof(real_T) << 4U);
    ipiv[0] = 1;
    ipiv[1] = 2;
    ipiv[2] = 3;
    ipiv[3] = 4;
    for (j = 0; j < 3; j++) {
        jBcol = j * 5 + 2;
        kAcol = j * 5;
        iy = 1;
        ix = jBcol - 2;
        temp = std::abs(A[kAcol]);
        for (k = 2; k <= 4 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > temp) {
                iy = k;
                temp = s;
            }
        }

        if (A[(jBcol + iy) - 3] != 0.0) {
            if (iy - 1 != 0) {
                iy += j;
                ipiv[j] = (int8_T) iy;
                iy--;
                temp = A[j];
                A[j] = A[iy];
                A[iy] = temp;
                ix = j + 4;
                iy += 4;
                temp = A[ix];
                A[ix] = A[iy];
                A[iy] = temp;
                ix += 4;
                iy += 4;
                temp = A[ix];
                A[ix] = A[iy];
                A[iy] = temp;
                ix += 4;
                iy += 4;
                temp = A[ix];
                A[ix] = A[iy];
                A[iy] = temp;
            }

            k = jBcol - j;
            for (ix = jBcol - 1; ix < k + 2; ix++) {
                A[ix] /= A[kAcol];
            }
        }

        iy = kAcol;
        kAcol += 4;
        for (j_0 = 1; j_0 <= 3 - j; j_0++) {
            if (A[kAcol] != 0.0) {
                temp = -A[kAcol];
                ix = jBcol - 1;
                k = iy - j;
                for (ijA = 5 + iy; ijA < k + 8; ijA++) {
                    A[ijA] += A[ix] * temp;
                    ix++;
                }
            }

            kAcol += 4;
            iy += 4;
        }
    }

    memcpy(&y[0], &u1[0], sizeof(real_T) << 4U);
    for (j = 0; j < 3; j++) {
        if (j + 1 != ipiv[j]) {
            jBcol = ipiv[j] - 1;
            temp = y[j];
            y[j] = y[jBcol];
            y[jBcol] = temp;
            temp = y[j + 4];
            y[j + 4] = y[jBcol + 4];
            y[jBcol + 4] = temp;
            temp = y[j + 8];
            y[j + 8] = y[jBcol + 8];
            y[jBcol + 8] = temp;
            temp = y[j + 12];
            y[j + 12] = y[jBcol + 12];
            y[jBcol + 12] = temp;
        }
    }

    for (j = 0; j < 4; j++) {
        jBcol = j << 2;
        if (y[jBcol] != 0.0) {
            for (ix = 1; ix + 1 < 5; ix++) {
                y[ix + jBcol] -= y[jBcol] * A[ix];
            }
        }

        if (y[1 + jBcol] != 0.0) {
            for (ix = 2; ix + 1 < 5; ix++) {
                y[ix + jBcol] -= y[1 + jBcol] * A[ix + 4];
            }
        }

        if (y[2 + jBcol] != 0.0) {
            for (ix = 3; ix + 1 < 5; ix++) {
                y[ix + jBcol] -= y[2 + jBcol] * A[ix + 8];
            }
        }
    }

    for (j = 0; j < 4; j++) {
        jBcol = j << 2;
        if (y[3 + jBcol] != 0.0) {
            y[3 + jBcol] /= A[15];
            for (ix = 0; ix < 3; ix++) {
                y[ix + jBcol] -= y[3 + jBcol] * A[ix + 12];
            }
        }

        if (y[2 + jBcol] != 0.0) {
            y[2 + jBcol] /= A[10];
            for (ix = 0; ix < 2; ix++) {
                y[ix + jBcol] -= y[2 + jBcol] * A[ix + 8];
            }
        }

        if (y[1 + jBcol] != 0.0) {
            y[1 + jBcol] /= A[5];
            for (ix = 0; ix < 1; ix++) {
                y[ix + jBcol] -= y[1 + jBcol] * A[ix + 4];
            }
        }

        if (y[jBcol] != 0.0) {
            y[jBcol] /= A[0];
        }
    }
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_eye(real_T I[16])
{
    memset(&I[0], 0, sizeof(real_T) << 4U);
    I[0] = 1.0;
    I[5] = 1.0;
    I[10] = 1.0;
    I[15] = 1.0;
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xgetrf_c(real_T A[16],
                                                                 int32_T ipiv[4],
                                                                 int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    ipiv[0] = 1;
    ipiv[1] = 2;
    ipiv[2] = 3;
    ipiv[3] = 4;
    *info = 0;
    for (j = 0; j < 3; j++) {
        c = j * 5;
        iy = 0;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 4 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                iy = k - 1;
                smax = s;
            }
        }

        if (A[c + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = iy + 1;
                smax = A[j];
                A[j] = A[iy];
                A[iy] = smax;
                ix = j + 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
                ix += 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
                ix += 4;
                iy += 4;
                smax = A[ix];
                A[ix] = A[iy];
                A[iy] = smax;
            }

            iy = (c - j) + 4;
            for (ix = c + 1; ix < iy; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        iy = c;
        ix = c + 4;
        for (k = 1; k <= 3 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                c_ix = c + 1;
                d = (iy - j) + 8;
                for (ijA = 5 + iy; ijA < d; ijA++) {
                    A[ijA] += A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 4;
            iy += 4;
        }
    }

    if ((*info == 0) && (!(A[15] != 0.0))) {
        *info = 4;
    }
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mrdivide_o(real_T A[16], const real_T B[16])
{
    real_T b_A[16];
    int32_T ipiv[4];
    int32_T info;
    real_T b_temp;
    int32_T jBcol;
    int32_T jAcol;
    int32_T kBcol;
    int32_T k;
    memcpy(&b_A[0], &B[0], sizeof(real_T) << 4U);
    torqueBalancingYoga_xgetrf_c(b_A, ipiv, &info);
    for (info = 0; info < 4; info++) {
        jBcol = info << 2;
        jAcol = info << 2;
        for (k = 1; k <= info; k++) {
            kBcol = (k - 1) << 2;
            if (b_A[(k + jAcol) - 1] != 0.0) {
                A[jBcol] -= b_A[(k + jAcol) - 1] * A[kBcol];
                A[1 + jBcol] -= b_A[(k + jAcol) - 1] * A[1 + kBcol];
                A[2 + jBcol] -= b_A[(k + jAcol) - 1] * A[2 + kBcol];
                A[3 + jBcol] -= b_A[(k + jAcol) - 1] * A[3 + kBcol];
            }
        }

        b_temp = 1.0 / b_A[info + jAcol];
        A[jBcol] *= b_temp;
        A[1 + jBcol] *= b_temp;
        A[2 + jBcol] *= b_temp;
        A[3 + jBcol] *= b_temp;
    }

    for (info = 3; info >= 0; info--) {
        jBcol = info << 2;
        jAcol = (info << 2) - 1;
        for (k = info + 2; k < 5; k++) {
            kBcol = (k - 1) << 2;
            if (b_A[k + jAcol] != 0.0) {
                A[jBcol] -= b_A[k + jAcol] * A[kBcol];
                A[1 + jBcol] -= b_A[k + jAcol] * A[1 + kBcol];
                A[2 + jBcol] -= b_A[k + jAcol] * A[2 + kBcol];
                A[3 + jBcol] -= b_A[k + jAcol] * A[3 + kBcol];
            }
        }
    }

    for (info = 2; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            jBcol = ipiv[info] - 1;
            b_temp = A[info << 2];
            A[info << 2] = A[jBcol << 2];
            A[jBcol << 2] = b_temp;
            b_temp = A[(info << 2) + 1];
            A[1 + (info << 2)] = A[(jBcol << 2) + 1];
            A[1 + (jBcol << 2)] = b_temp;
            b_temp = A[(info << 2) + 2];
            A[2 + (info << 2)] = A[(jBcol << 2) + 2];
            A[2 + (jBcol << 2)] = b_temp;
            b_temp = A[(info << 2) + 3];
            A[3 + (info << 2)] = A[(jBcol << 2) + 3];
            A[3 + (jBcol << 2)] = b_temp;
        }
    }
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mldivide_e(const real_T A[16], real_T B[4])
{
    real_T temp;
    real_T b_A[16];
    int32_T ipiv[4];
    int32_T info;
    memcpy(&b_A[0], &A[0], sizeof(real_T) << 4U);
    torqueBalancingYoga_xgetrf_c(b_A, ipiv, &info);
    if (ipiv[0] != 1) {
        temp = B[0];
        B[0] = B[ipiv[0] - 1];
        B[ipiv[0] - 1] = temp;
    }

    if (ipiv[1] != 2) {
        temp = B[1];
        B[1] = B[ipiv[1] - 1];
        B[ipiv[1] - 1] = temp;
    }

    if (ipiv[2] != 3) {
        temp = B[2];
        B[2] = B[ipiv[2] - 1];
        B[ipiv[2] - 1] = temp;
    }

    if (B[0] != 0.0) {
        for (info = 1; info + 1 < 5; info++) {
            B[info] -= B[0] * b_A[info];
        }
    }

    if (B[1] != 0.0) {
        for (info = 2; info + 1 < 5; info++) {
            B[info] -= b_A[info + 4] * B[1];
        }
    }

    if (B[2] != 0.0) {
        for (info = 3; info + 1 < 5; info++) {
            B[info] -= b_A[info + 8] * B[2];
        }
    }

    if (B[3] != 0.0) {
        B[3] /= b_A[15];
        for (info = 0; info < 3; info++) {
            B[info] -= b_A[info + 12] * B[3];
        }
    }

    if (B[2] != 0.0) {
        B[2] /= b_A[10];
        for (info = 0; info < 2; info++) {
            B[info] -= b_A[info + 8] * B[2];
        }
    }

    if (B[1] != 0.0) {
        B[1] /= b_A[5];
        for (info = 0; info < 1; info++) {
            B[info] -= b_A[info + 4] * B[1];
        }
    }

    if (B[0] != 0.0) {
        B[0] /= b_A[0];
    }
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_norm(const real_T x[6])
{
    real_T y;
    real_T scale;
    real_T absxk;
    real_T t;
    int32_T k;
    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (k = 0; k < 6; k++) {
        absxk = std::abs(x[k]);
        if (absxk > scale) {
            t = scale / absxk;
            y = y * t * t + 1.0;
            scale = absxk;
        }
        else {
            t = absxk / scale;
            y += t * t;
        }
    }

    return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S34>/stateMachineYogaFCN' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_norm_i(const real_T x[2])
{
    real_T y;
    real_T scale;
    real_T absxk;
    real_T t;
    scale = 3.3121686421112381E-170;
    absxk = std::abs(x[0]);
    if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
    }
    else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
    }

    absxk = std::abs(x[1]);
    if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
    }
    else {
        t = absxk / scale;
        y += t * t;
    }

    return scale * std::sqrt(y);
}

/* Function for MATLAB Function: '<S36>/Compute Base Velocity' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mldivide(const real_T A[36], real_T B[72])
{
    int32_T ip;
    real_T b_A[36];
    int8_T ipiv[6];
    int32_T j;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    memcpy(&b_A[0], &A[0], 36U * sizeof(real_T));
    for (c_ix = 0; c_ix < 6; c_ix++) {
        ipiv[c_ix] = (int8_T)(1 + c_ix);
    }

    for (j = 0; j < 5; j++) {
        ip = j * 7;
        iy = 0;
        ix = ip;
        smax = std::abs(b_A[ip]);
        for (k = 2; k <= 6 - j; k++) {
            ix++;
            s = std::abs(b_A[ix]);
            if (s > smax) {
                iy = k - 1;
                smax = s;
            }
        }

        if (b_A[ip + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = (int8_T)(iy + 1);
                ix = j;
                for (k = 0; k < 6; k++) {
                    smax = b_A[ix];
                    b_A[ix] = b_A[iy];
                    b_A[iy] = smax;
                    ix += 6;
                    iy += 6;
                }
            }

            iy = (ip - j) + 6;
            for (ix = ip + 1; ix < iy; ix++) {
                b_A[ix] /= b_A[ip];
            }
        }

        iy = ip;
        ix = ip + 6;
        for (k = 1; k <= 5 - j; k++) {
            smax = b_A[ix];
            if (b_A[ix] != 0.0) {
                c_ix = ip + 1;
                d = (iy - j) + 12;
                for (ijA = 7 + iy; ijA < d; ijA++) {
                    b_A[ijA] += b_A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 6;
            iy += 6;
        }
    }

    for (j = 0; j < 5; j++) {
        if (j + 1 != ipiv[j]) {
            ip = ipiv[j] - 1;
            for (iy = 0; iy < 12; iy++) {
                smax = B[6 * iy + j];
                B[j + 6 * iy] = B[6 * iy + ip];
                B[ip + 6 * iy] = smax;
            }
        }
    }

    for (j = 0; j < 12; j++) {
        ip = 6 * j;
        for (iy = 0; iy < 6; iy++) {
            ix = 6 * iy;
            if (B[iy + ip] != 0.0) {
                for (k = iy + 1; k + 1 < 7; k++) {
                    c_ix = k + ip;
                    B[c_ix] -= B[iy + ip] * b_A[k + ix];
                }
            }
        }
    }

    for (j = 0; j < 12; j++) {
        ip = 6 * j;
        for (iy = 5; iy >= 0; iy--) {
            ix = 6 * iy;
            if (B[iy + ip] != 0.0) {
                B[iy + ip] /= b_A[iy + ix];
                for (k = 0; k < iy; k++) {
                    c_ix = k + ip;
                    B[c_ix] -= B[iy + ip] * b_A[k + ix];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S15>/Add motor reflected inertias' */
void torqueBalancingYogaModelClass::torqueBala_computeMotorsInertia(
    const struct_JFJmPMlc0Zv9CJd3pEkJyG* b_Config,
    real_T reflectedInertia[529])
{
    static const real_T b[529] = {
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

    int32_T i;
    int32_T i_0;
    int32_T i_1;
    int32_T b_Config_tmp;
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (i = 0; i < 23; i++) {
            b_Config_tmp = i + 23 * i_0;
            torqueBalancingYoga_B.b_Config[b_Config_tmp] = 0.0;
            for (i_1 = 0; i_1 < 23; i_1++) {
                torqueBalancingYoga_B.b_Config[b_Config_tmp] =
                    b_Config->T[23 * i_1 + i] * b_Config->Gamma[23 * i_0 + i_1]
                    + torqueBalancingYoga_B.b_Config[23 * i_0 + i];
            }
        }
    }

    torqueBalancingYoga_mrdivide(
        b, torqueBalancingYoga_B.b_Config, torqueBalancingYoga_B.invTGamma_c);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (i = 0; i < 23; i++) {
            b_Config_tmp = i_0 + 23 * i;
            torqueBalancingYoga_B.b_Config[b_Config_tmp] = 0.0;
            for (i_1 = 0; i_1 < 23; i_1++) {
                torqueBalancingYoga_B.b_Config[b_Config_tmp] =
                    b_Config->T[23 * i_1 + i] * b_Config->Gamma[23 * i_0 + i_1]
                    + torqueBalancingYoga_B.b_Config[23 * i + i_0];
            }
        }
    }

    torqueBalancingYoga_mrdivide(
        b, torqueBalancingYoga_B.b_Config, torqueBalancingYoga_B.invTGamma_t);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (i = 0; i < 23; i++) {
            b_Config_tmp = i_0 + 23 * i;
            torqueBalancingYoga_B.b_Config[b_Config_tmp] = 0.0;
            for (i_1 = 0; i_1 < 23; i_1++) {
                torqueBalancingYoga_B.b_Config[b_Config_tmp] =
                    torqueBalancingYoga_B.invTGamma_t[23 * i_1 + i_0] * b_Config->I_m[23 * i + i_1]
                    + torqueBalancingYoga_B.b_Config[23 * i + i_0];
            }
        }

        for (i = 0; i < 23; i++) {
            b_Config_tmp = i_0 + 23 * i;
            reflectedInertia[b_Config_tmp] = 0.0;
            for (i_1 = 0; i_1 < 23; i_1++) {
                reflectedInertia[b_Config_tmp] =
                    torqueBalancingYoga_B.b_Config[23 * i_1 + i_0]
                        * torqueBalancingYoga_B.invTGamma_c[23 * i + i_1]
                    + reflectedInertia[23 * i + i_0];
            }
        }
    }
}

/* Function for MATLAB Function: '<S109>/References for L' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mldivide_h(const real_T A[36], real_T B[36])
{
    int32_T ip;
    real_T b_A[36];
    int8_T ipiv[6];
    int32_T j;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    memcpy(&b_A[0], &A[0], 36U * sizeof(real_T));
    for (c_ix = 0; c_ix < 6; c_ix++) {
        ipiv[c_ix] = (int8_T)(1 + c_ix);
    }

    for (j = 0; j < 5; j++) {
        ip = j * 7;
        iy = 0;
        ix = ip;
        smax = std::abs(b_A[ip]);
        for (k = 2; k <= 6 - j; k++) {
            ix++;
            s = std::abs(b_A[ix]);
            if (s > smax) {
                iy = k - 1;
                smax = s;
            }
        }

        if (b_A[ip + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = (int8_T)(iy + 1);
                ix = j;
                for (k = 0; k < 6; k++) {
                    smax = b_A[ix];
                    b_A[ix] = b_A[iy];
                    b_A[iy] = smax;
                    ix += 6;
                    iy += 6;
                }
            }

            iy = (ip - j) + 6;
            for (ix = ip + 1; ix < iy; ix++) {
                b_A[ix] /= b_A[ip];
            }
        }

        iy = ip;
        ix = ip + 6;
        for (k = 1; k <= 5 - j; k++) {
            smax = b_A[ix];
            if (b_A[ix] != 0.0) {
                c_ix = ip + 1;
                d = (iy - j) + 12;
                for (ijA = 7 + iy; ijA < d; ijA++) {
                    b_A[ijA] += b_A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 6;
            iy += 6;
        }
    }

    for (j = 0; j < 5; j++) {
        if (j + 1 != ipiv[j]) {
            ip = ipiv[j] - 1;
            for (iy = 0; iy < 6; iy++) {
                smax = B[6 * iy + j];
                B[j + 6 * iy] = B[6 * iy + ip];
                B[ip + 6 * iy] = smax;
            }
        }
    }

    for (j = 0; j < 6; j++) {
        ip = 6 * j;
        for (iy = 0; iy < 6; iy++) {
            ix = 6 * iy;
            if (B[iy + ip] != 0.0) {
                for (k = iy + 1; k + 1 < 7; k++) {
                    c_ix = k + ip;
                    B[c_ix] -= B[iy + ip] * b_A[k + ix];
                }
            }
        }
    }

    for (j = 0; j < 6; j++) {
        ip = 6 * j;
        for (iy = 5; iy >= 0; iy--) {
            ix = 6 * iy;
            if (B[iy + ip] != 0.0) {
                B[iy + ip] /= b_A[iy + ix];
                for (k = 0; k < iy; k++) {
                    c_ix = k + ip;
                    B[c_ix] -= B[iy + ip] * b_A[k + ix];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xzgetrf(real_T A[36],
                                                                int32_T ipiv[6],
                                                                int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    for (j = 0; j < 6; j++) {
        ipiv[j] = 1 + j;
    }

    *info = 0;
    for (j = 0; j < 5; j++) {
        c = j * 7;
        iy = 0;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 6 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                iy = k - 1;
                smax = s;
            }
        }

        if (A[c + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = iy + 1;
                ix = j;
                for (k = 0; k < 6; k++) {
                    smax = A[ix];
                    A[ix] = A[iy];
                    A[iy] = smax;
                    ix += 6;
                    iy += 6;
                }
            }

            iy = (c - j) + 6;
            for (ix = c + 1; ix < iy; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        iy = c;
        ix = c + 6;
        for (k = 1; k <= 5 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                c_ix = c + 1;
                d = (iy - j) + 12;
                for (ijA = 7 + iy; ijA < d; ijA++) {
                    A[ijA] += A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 6;
            iy += 6;
        }
    }

    if ((*info == 0) && (!(A[35] != 0.0))) {
        *info = 6;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_inv(const real_T x[36], real_T y[36])
{
    int8_T p[6];
    int32_T j;
    real_T A[36];
    int32_T ipiv[6];
    int32_T info;
    int32_T pipk;
    int32_T kAcol;
    int32_T b_i;
    int32_T y_tmp;
    for (b_i = 0; b_i < 36; b_i++) {
        y[b_i] = 0.0;
        A[b_i] = x[b_i];
    }

    torqueBalancingYoga_xzgetrf(A, ipiv, &info);
    for (b_i = 0; b_i < 6; b_i++) {
        p[b_i] = (int8_T)(1 + b_i);
    }

    for (info = 0; info < 5; info++) {
        if (ipiv[info] > 1 + info) {
            pipk = p[ipiv[info] - 1];
            p[ipiv[info] - 1] = p[info];
            p[info] = (int8_T) pipk;
        }
    }

    for (info = 0; info < 6; info++) {
        pipk = p[info] - 1;
        y[info + 6 * pipk] = 1.0;
        for (j = info; j + 1 < 7; j++) {
            if (y[6 * pipk + j] != 0.0) {
                for (kAcol = j + 1; kAcol + 1 < 7; kAcol++) {
                    y[kAcol + 6 * pipk] -= y[6 * pipk + j] * A[6 * j + kAcol];
                }
            }
        }
    }

    for (info = 0; info < 6; info++) {
        pipk = 6 * info;
        for (j = 5; j >= 0; j--) {
            kAcol = 6 * j;
            b_i = j + pipk;
            if (y[b_i] != 0.0) {
                y[b_i] = y[j + pipk] / A[j + kAcol];
                for (b_i = 0; b_i < j; b_i++) {
                    y_tmp = b_i + pipk;
                    y[y_tmp] -= y[j + pipk] * A[b_i + kAcol];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T
torqueBalancingYogaModelClass::torqueBalancingYoga_xnrm2(int32_T n, const real_T x[72], int32_T ix0)
{
    real_T y;
    real_T scale;
    int32_T kend;
    real_T absxk;
    real_T t;
    int32_T k;
    y = 0.0;
    if (!(n < 1)) {
        if (n == 1) {
            y = std::abs(x[ix0 - 1]);
        }
        else {
            scale = 3.3121686421112381E-170;
            kend = (ix0 + n) - 1;
            for (k = ix0; k <= kend; k++) {
                absxk = std::abs(x[k - 1]);
                if (absxk > scale) {
                    t = scale / absxk;
                    y = y * t * t + 1.0;
                    scale = absxk;
                }
                else {
                    t = absxk / scale;
                    y += t * t;
                }
            }

            y = scale * std::sqrt(y);
        }
    }

    return y;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xnrm2_m(int32_T n,
                                                                  const real_T x[6],
                                                                  int32_T ix0)
{
    real_T y;
    real_T scale;
    int32_T kend;
    real_T absxk;
    real_T t;
    int32_T k;
    y = 0.0;
    if (!(n < 1)) {
        if (n == 1) {
            y = std::abs(x[ix0 - 1]);
        }
        else {
            scale = 3.3121686421112381E-170;
            kend = (ix0 + n) - 1;
            for (k = ix0; k <= kend; k++) {
                absxk = std::abs(x[k - 1]);
                if (absxk > scale) {
                    t = scale / absxk;
                    y = y * t * t + 1.0;
                    scale = absxk;
                }
                else {
                    t = absxk / scale;
                    y += t * t;
                }
            }

            y = scale * std::sqrt(y);
        }
    }

    return y;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy_er(int32_T n,
                                                                 real_T a,
                                                                 const real_T x[12],
                                                                 int32_T ix0,
                                                                 real_T y[72],
                                                                 int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * x[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy_e(int32_T n,
                                                                real_T a,
                                                                const real_T x[72],
                                                                int32_T ix0,
                                                                real_T y[12],
                                                                int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * x[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xdotc(int32_T n,
                                                                const real_T x[72],
                                                                int32_T ix0,
                                                                const real_T y[72],
                                                                int32_T iy0)
{
    real_T d;
    int32_T ix;
    int32_T iy;
    int32_T k;
    d = 0.0;
    if (!(n < 1)) {
        ix = ix0;
        iy = iy0;
        for (k = 1; k <= n; k++) {
            d += x[ix - 1] * y[iy - 1];
            ix++;
            iy++;
        }
    }

    return d;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy(int32_T n,
                                                              real_T a,
                                                              int32_T ix0,
                                                              real_T y[72],
                                                              int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * y[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xdotc_l(int32_T n,
                                                                  const real_T x[36],
                                                                  int32_T ix0,
                                                                  const real_T y[36],
                                                                  int32_T iy0)
{
    real_T d;
    int32_T ix;
    int32_T iy;
    int32_T k;
    d = 0.0;
    if (!(n < 1)) {
        ix = ix0;
        iy = iy0;
        for (k = 1; k <= n; k++) {
            d += x[ix - 1] * y[iy - 1];
            ix++;
            iy++;
        }
    }

    return d;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy_erc(int32_T n,
                                                                  real_T a,
                                                                  int32_T ix0,
                                                                  real_T y[36],
                                                                  int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * y[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xscal(real_T a, real_T x[72], int32_T ix0)
{
    int32_T k;
    for (k = ix0; k <= ix0 + 11; k++) {
        x[k - 1] *= a;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xscal_m(real_T a, real_T x[36], int32_T ix0)
{
    int32_T k;
    for (k = ix0; k <= ix0 + 5; k++) {
        x[k - 1] *= a;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xswap_m(real_T x[36],
                                                                int32_T ix0,
                                                                int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 6; k++) {
        temp = x[ix];
        x[ix] = x[iy];
        x[iy] = temp;
        ix++;
        iy++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xswap_mc(real_T x[72],
                                                                 int32_T ix0,
                                                                 int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 12; k++) {
        temp = x[ix];
        x[ix] = x[iy];
        x[iy] = temp;
        ix++;
        iy++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xrotg(real_T* a,
                                                              real_T* b,
                                                              real_T* c,
                                                              real_T* s)
{
    real_T roe;
    real_T absa;
    real_T absb;
    real_T scale;
    real_T ads;
    real_T bds;
    roe = *b;
    absa = std::abs(*a);
    absb = std::abs(*b);
    if (absa > absb) {
        roe = *a;
    }

    scale = absa + absb;
    if (scale == 0.0) {
        *s = 0.0;
        *c = 1.0;
        absa = 0.0;
    }
    else {
        ads = absa / scale;
        bds = absb / scale;
        scale *= std::sqrt(ads * ads + bds * bds);
        if (roe < 0.0) {
            scale = -scale;
        }

        *c = *a / scale;
        *s = *b / scale;
        if (absa > absb) {
            absa = *s;
        }
        else if (*c != 0.0) {
            absa = 1.0 / *c;
        }
        else {
            absa = 1.0;
        }
    }

    *a = scale;
    *b = absa;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xrot(real_T x[36],
                                                             int32_T ix0,
                                                             int32_T iy0,
                                                             real_T c,
                                                             real_T s)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 6; k++) {
        temp = c * x[ix] + s * x[iy];
        x[iy] = c * x[iy] - s * x[ix];
        x[ix] = temp;
        iy++;
        ix++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xrot_g(real_T x[72],
                                                               int32_T ix0,
                                                               int32_T iy0,
                                                               real_T c,
                                                               real_T s)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 12; k++) {
        temp = c * x[ix] + s * x[iy];
        x[iy] = c * x[iy] - s * x[ix];
        x[ix] = temp;
        iy++;
        ix++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_svd(const real_T A[72],
                                                            real_T U[72],
                                                            real_T s[6],
                                                            real_T V[36])
{
    real_T b_A[72];
    real_T b_s[6];
    real_T e[6];
    real_T work[12];
    real_T Vf[36];
    int32_T qq;
    boolean_T apply_transform;
    real_T nrm;
    int32_T qjj;
    int32_T qp1jj;
    int32_T qp1q;
    real_T rt;
    real_T ztest;
    real_T smm1;
    real_T emm1;
    real_T sqds;
    real_T shift;
    int32_T j_ii;
    int32_T i;
    int32_T exitg1;
    boolean_T exitg2;
    memcpy(&b_A[0], &A[0], 72U * sizeof(real_T));
    for (i = 0; i < 6; i++) {
        b_s[i] = 0.0;
        e[i] = 0.0;
    }

    memset(&work[0], 0, 12U * sizeof(real_T));
    memset(&U[0], 0, 72U * sizeof(real_T));
    memset(&Vf[0], 0, 36U * sizeof(real_T));
    for (i = 0; i < 6; i++) {
        qq = 12 * i + i;
        apply_transform = false;
        nrm = torqueBalancingYoga_xnrm2(12 - i, b_A, qq + 1);
        if (nrm > 0.0) {
            apply_transform = true;
            if (b_A[qq] < 0.0) {
                b_s[i] = -nrm;
            }
            else {
                b_s[i] = nrm;
            }

            if (std::abs(b_s[i]) >= 1.0020841800044864E-292) {
                nrm = 1.0 / b_s[i];
                qp1q = (qq - i) + 12;
                for (qjj = qq; qjj < qp1q; qjj++) {
                    b_A[qjj] *= nrm;
                }
            }
            else {
                qp1q = (qq - i) + 12;
                for (qjj = qq; qjj < qp1q; qjj++) {
                    b_A[qjj] /= b_s[i];
                }
            }

            b_A[qq]++;
            b_s[i] = -b_s[i];
        }
        else {
            b_s[i] = 0.0;
        }

        for (qp1q = i + 1; qp1q + 1 < 7; qp1q++) {
            qjj = 12 * qp1q + i;
            if (apply_transform) {
                torqueBalancingYoga_xaxpy(
                    12 - i,
                    -(torqueBalancingYoga_xdotc(12 - i, b_A, qq + 1, b_A, qjj + 1)
                      / b_A[i + 12 * i]),
                    qq + 1,
                    b_A,
                    qjj + 1);
            }

            e[qp1q] = b_A[qjj];
        }

        for (qq = i; qq + 1 < 13; qq++) {
            U[qq + 12 * i] = b_A[12 * i + qq];
        }

        if (i + 1 <= 4) {
            nrm = torqueBalancingYoga_xnrm2_m(5 - i, e, i + 2);
            if (nrm == 0.0) {
                e[i] = 0.0;
            }
            else {
                if (e[i + 1] < 0.0) {
                    e[i] = -nrm;
                }
                else {
                    e[i] = nrm;
                }

                nrm = e[i];
                if (std::abs(e[i]) >= 1.0020841800044864E-292) {
                    nrm = 1.0 / e[i];
                    for (qq = i + 1; qq < 6; qq++) {
                        e[qq] *= nrm;
                    }
                }
                else {
                    for (qq = i + 1; qq < 6; qq++) {
                        e[qq] /= nrm;
                    }
                }

                e[i + 1]++;
                e[i] = -e[i];
                for (qq = i + 1; qq + 1 < 13; qq++) {
                    work[qq] = 0.0;
                }

                for (qq = i + 1; qq + 1 < 7; qq++) {
                    torqueBalancingYoga_xaxpy_e(11 - i, e[qq], b_A, (i + 12 * qq) + 2, work, i + 2);
                }

                for (qq = i + 1; qq + 1 < 7; qq++) {
                    torqueBalancingYoga_xaxpy_er(
                        11 - i, -e[qq] / e[i + 1], work, i + 2, b_A, (i + 12 * qq) + 2);
                }
            }

            for (qq = i + 1; qq + 1 < 7; qq++) {
                Vf[qq + 6 * i] = e[qq];
            }
        }
    }

    i = 4;
    e[4] = b_A[64];
    e[5] = 0.0;
    for (qp1q = 5; qp1q >= 0; qp1q--) {
        qq = 12 * qp1q + qp1q;
        if (b_s[qp1q] != 0.0) {
            for (qp1jj = qp1q + 1; qp1jj + 1 < 7; qp1jj++) {
                qjj = (12 * qp1jj + qp1q) + 1;
                torqueBalancingYoga_xaxpy(
                    12 - qp1q,
                    -(torqueBalancingYoga_xdotc(12 - qp1q, U, qq + 1, U, qjj) / U[qq]),
                    qq + 1,
                    U,
                    qjj);
            }

            for (qjj = qp1q; qjj + 1 < 13; qjj++) {
                U[qjj + 12 * qp1q] = -U[12 * qp1q + qjj];
            }

            U[qq]++;
            for (qq = 1; qq <= qp1q; qq++) {
                U[(qq + 12 * qp1q) - 1] = 0.0;
            }
        }
        else {
            memset(&U[qp1q * 12], 0, 12U * sizeof(real_T));
            U[qq] = 1.0;
        }
    }

    for (qq = 5; qq >= 0; qq--) {
        if ((qq + 1 <= 4) && (e[qq] != 0.0)) {
            qp1q = (6 * qq + qq) + 2;
            for (qjj = qq + 1; qjj + 1 < 7; qjj++) {
                qp1jj = (6 * qjj + qq) + 2;
                torqueBalancingYoga_xaxpy_erc(
                    5 - qq,
                    -(torqueBalancingYoga_xdotc_l(5 - qq, Vf, qp1q, Vf, qp1jj) / Vf[qp1q - 1]),
                    qp1q,
                    Vf,
                    qp1jj);
            }
        }

        for (qp1q = 0; qp1q < 6; qp1q++) {
            Vf[qp1q + 6 * qq] = 0.0;
        }

        Vf[qq + 6 * qq] = 1.0;
    }

    for (qq = 0; qq < 6; qq++) {
        ztest = e[qq];
        if (b_s[qq] != 0.0) {
            rt = std::abs(b_s[qq]);
            nrm = b_s[qq] / rt;
            b_s[qq] = rt;
            if (qq + 1 < 6) {
                ztest = e[qq] / nrm;
            }

            torqueBalancingYoga_xscal(nrm, U, 1 + 12 * qq);
        }

        if ((qq + 1 < 6) && (ztest != 0.0)) {
            rt = std::abs(ztest);
            nrm = rt / ztest;
            ztest = rt;
            b_s[qq + 1] *= nrm;
            torqueBalancingYoga_xscal_m(nrm, Vf, 1 + 6 * (qq + 1));
        }

        e[qq] = ztest;
    }

    qq = 0;
    nrm = 0.0;
    for (qp1q = 0; qp1q < 6; qp1q++) {
        ztest = std::abs(b_s[qp1q]);
        rt = std::abs(e[qp1q]);
        if (ztest > rt) {
            rt = ztest;
        }

        if (!(nrm > rt)) {
            nrm = rt;
        }
    }

    while ((i + 2 > 0) && (!(qq >= 75))) {
        qp1jj = i + 1;
        do {
            exitg1 = 0;
            qp1q = qp1jj;
            if (qp1jj == 0) {
                exitg1 = 1;
            }
            else {
                rt = std::abs(e[qp1jj - 1]);
                if ((rt
                     <= (std::abs(b_s[qp1jj - 1]) + std::abs(b_s[qp1jj])) * 2.2204460492503131E-16)
                    || ((rt <= 1.0020841800044864E-292)
                        || ((qq > 20) && (rt <= 2.2204460492503131E-16 * nrm)))) {
                    e[qp1jj - 1] = 0.0;
                    exitg1 = 1;
                }
                else {
                    qp1jj--;
                }
            }
        } while (exitg1 == 0);

        if (i + 1 == qp1jj) {
            qp1jj = 4;
        }
        else {
            qjj = i + 2;
            j_ii = i + 2;
            exitg2 = false;
            while ((!exitg2) && (j_ii >= qp1jj)) {
                qjj = j_ii;
                if (j_ii == qp1jj) {
                    exitg2 = true;
                }
                else {
                    rt = 0.0;
                    if (j_ii < i + 2) {
                        rt = std::abs(e[j_ii - 1]);
                    }

                    if (j_ii > qp1jj + 1) {
                        rt += std::abs(e[j_ii - 2]);
                    }

                    ztest = std::abs(b_s[j_ii - 1]);
                    if ((ztest <= 2.2204460492503131E-16 * rt)
                        || (ztest <= 1.0020841800044864E-292)) {
                        b_s[j_ii - 1] = 0.0;
                        exitg2 = true;
                    }
                    else {
                        j_ii--;
                    }
                }
            }

            if (qjj == qp1jj) {
                qp1jj = 3;
            }
            else if (i + 2 == qjj) {
                qp1jj = 1;
            }
            else {
                qp1jj = 2;
                qp1q = qjj;
            }
        }

        switch (qp1jj) {
            case 1:
                rt = e[i];
                e[i] = 0.0;
                for (qjj = i; qjj + 1 >= qp1q + 1; qjj--) {
                    ztest = b_s[qjj];
                    torqueBalancingYoga_xrotg(&ztest, &rt, &sqds, &smm1);
                    b_s[qjj] = ztest;
                    if (qjj + 1 > qp1q + 1) {
                        rt = e[qjj - 1] * -smm1;
                        e[qjj - 1] *= sqds;
                    }

                    torqueBalancingYoga_xrot(Vf, 1 + 6 * qjj, 1 + 6 * (i + 1), sqds, smm1);
                }
                break;

            case 2:
                rt = e[qp1q - 1];
                e[qp1q - 1] = 0.0;
                for (qjj = qp1q; qjj < i + 2; qjj++) {
                    ztest = b_s[qjj];
                    torqueBalancingYoga_xrotg(&ztest, &rt, &sqds, &smm1);
                    b_s[qjj] = ztest;
                    rt = -smm1 * e[qjj];
                    e[qjj] *= sqds;
                    torqueBalancingYoga_xrot_g(U, 1 + 12 * qjj, 1 + 12 * (qp1q - 1), sqds, smm1);
                }
                break;

            case 3:
                ztest = std::abs(b_s[i + 1]);
                rt = std::abs(b_s[i]);
                if (ztest > rt) {
                    rt = ztest;
                }

                ztest = std::abs(e[i]);
                if (rt > ztest) {
                    ztest = rt;
                }

                rt = std::abs(b_s[qp1q]);
                if (ztest > rt) {
                    rt = ztest;
                }

                ztest = std::abs(e[qp1q]);
                if (rt > ztest) {
                    ztest = rt;
                }

                rt = b_s[i + 1] / ztest;
                smm1 = b_s[i] / ztest;
                emm1 = e[i] / ztest;
                sqds = b_s[qp1q] / ztest;
                smm1 = ((smm1 + rt) * (smm1 - rt) + emm1 * emm1) / 2.0;
                emm1 *= rt;
                emm1 *= emm1;
                if ((smm1 != 0.0) || (emm1 != 0.0)) {
                    shift = std::sqrt(smm1 * smm1 + emm1);
                    if (smm1 < 0.0) {
                        shift = -shift;
                    }

                    shift = emm1 / (smm1 + shift);
                }
                else {
                    shift = 0.0;
                }

                rt = (sqds + rt) * (sqds - rt) + shift;
                ztest = e[qp1q] / ztest * sqds;
                for (qjj = qp1q + 1; qjj <= i + 1; qjj++) {
                    torqueBalancingYoga_xrotg(&rt, &ztest, &sqds, &smm1);
                    if (qjj > qp1q + 1) {
                        e[qjj - 2] = rt;
                    }

                    rt = b_s[qjj - 1] * sqds + e[qjj - 1] * smm1;
                    e[qjj - 1] = e[qjj - 1] * sqds - b_s[qjj - 1] * smm1;
                    ztest = smm1 * b_s[qjj];
                    b_s[qjj] *= sqds;
                    torqueBalancingYoga_xrot(Vf, 1 + 6 * (qjj - 1), 1 + 6 * qjj, sqds, smm1);
                    torqueBalancingYoga_xrotg(&rt, &ztest, &sqds, &smm1);
                    b_s[qjj - 1] = rt;
                    rt = e[qjj - 1] * sqds + smm1 * b_s[qjj];
                    b_s[qjj] = e[qjj - 1] * -smm1 + sqds * b_s[qjj];
                    ztest = smm1 * e[qjj];
                    e[qjj] *= sqds;
                    torqueBalancingYoga_xrot_g(U, 1 + 12 * (qjj - 1), 1 + 12 * qjj, sqds, smm1);
                }

                e[i] = rt;
                qq++;
                break;

            default:
                if (b_s[qp1q] < 0.0) {
                    b_s[qp1q] = -b_s[qp1q];
                    torqueBalancingYoga_xscal_m(-1.0, Vf, 1 + 6 * qp1q);
                }

                qq = qp1q + 1;
                while ((qp1q + 1 < 6) && (b_s[qp1q] < b_s[qq])) {
                    rt = b_s[qp1q];
                    b_s[qp1q] = b_s[qq];
                    b_s[qq] = rt;
                    torqueBalancingYoga_xswap_m(Vf, 1 + 6 * qp1q, 1 + 6 * (qp1q + 1));
                    torqueBalancingYoga_xswap_mc(U, 1 + 12 * qp1q, 1 + 12 * (qp1q + 1));
                    qp1q = qq;
                    qq++;
                }

                qq = 0;
                i--;
                break;
        }
    }

    for (i = 0; i < 6; i++) {
        s[i] = b_s[i];
        for (qq = 0; qq < 6; qq++) {
            V[qq + 6 * i] = Vf[6 * i + qq];
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_pinv(const real_T A[72],
                                                             real_T tol,
                                                             real_T X[72])
{
    real_T b_X[72];
    real_T V[36];
    int32_T r;
    int32_T vcol;
    real_T U[72];
    real_T s[6];
    int32_T j;
    int32_T ar;
    int32_T ia;
    int32_T b;
    int32_T ib;
    int32_T b_ic;
    real_T z;
    real_T A_0[72];
    memset(&b_X[0], 0, 72U * sizeof(real_T));
    for (r = 0; r < 6; r++) {
        for (vcol = 0; vcol < 12; vcol++) {
            A_0[vcol + 12 * r] = A[6 * vcol + r];
        }
    }

    torqueBalancingYoga_svd(A_0, U, s, V);
    r = 0;
    vcol = 1;
    while ((vcol < 7) && (s[vcol - 1] > tol)) {
        r++;
        vcol++;
    }

    if (r > 0) {
        vcol = 0;
        for (j = 1; j <= r; j++) {
            z = 1.0 / s[j - 1];
            for (ar = vcol; ar < vcol + 6; ar++) {
                V[ar] *= z;
            }

            vcol += 6;
        }

        for (vcol = 0; vcol <= 67; vcol += 6) {
            for (j = vcol; j < vcol + 6; j++) {
                b_X[j] = 0.0;
            }
        }

        vcol = -1;
        for (j = 0; j <= 67; j += 6) {
            ar = -1;
            vcol++;
            b = ((r - 1) * 12 + vcol) + 1;
            for (ib = vcol; ib + 1 <= b; ib += 12) {
                if (U[ib] != 0.0) {
                    ia = ar;
                    for (b_ic = j; b_ic < j + 6; b_ic++) {
                        ia++;
                        b_X[b_ic] += U[ib] * V[ia];
                    }
                }

                ar += 6;
            }
        }
    }

    for (r = 0; r < 6; r++) {
        for (vcol = 0; vcol < 12; vcol++) {
            X[vcol + 12 * r] = b_X[6 * vcol + r];
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_eye_a(real_T I[144])
{
    int32_T k;
    memset(&I[0], 0, 144U * sizeof(real_T));
    for (k = 0; k < 12; k++) {
        I[k + 12 * k] = 1.0;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xswap_mc0(real_T x[841],
                                                                  int32_T ix0,
                                                                  int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 29; k++) {
        temp = x[ix];
        x[ix] = x[iy];
        x[iy] = temp;
        ix += 29;
        iy += 29;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xgetrf_i(real_T A[841],
                                                                 int32_T ipiv[29],
                                                                 int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T jA;
    int32_T b_ix;
    int32_T d;
    int32_T ijA;
    for (j = 0; j < 29; j++) {
        ipiv[j] = 1 + j;
    }

    *info = 0;
    for (j = 0; j < 28; j++) {
        c = j * 30;
        jA = 1;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 29 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                jA = k;
                smax = s;
            }
        }

        if (A[(c + jA) - 1] != 0.0) {
            if (jA - 1 != 0) {
                ipiv[j] = j + jA;
                torqueBalancingYoga_xswap_mc0(A, j + 1, j + jA);
            }

            jA = (c - j) + 29;
            for (ix = c + 1; ix < jA; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        jA = c;
        ix = c + 29;
        for (k = 1; k <= 28 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                b_ix = c + 1;
                d = (jA - j) + 58;
                for (ijA = 30 + jA; ijA < d; ijA++) {
                    A[ijA] += A[b_ix] * -smax;
                    b_ix++;
                }
            }

            ix += 29;
            jA += 29;
        }
    }

    if ((*info == 0) && (!(A[840] != 0.0))) {
        *info = 29;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm_o(const real_T A[841], real_T B[348])
{
    real_T temp;
    int32_T jBcol;
    int32_T jAcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 0; j < 29; j++) {
        jBcol = 12 * j;
        jAcol = 29 * j;
        for (k = 1; k <= j; k++) {
            kBcol = (k - 1) * 12;
            if (A[(k + jAcol) - 1] != 0.0) {
                for (i = 0; i < 12; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[(k + jAcol) - 1] * B[i + kBcol];
                }
            }
        }

        temp = 1.0 / A[j + jAcol];
        for (jAcol = 0; jAcol < 12; jAcol++) {
            tmp = jAcol + jBcol;
            B[tmp] *= temp;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm_oz(const real_T A[841], real_T B[348])
{
    int32_T jAcol;
    int32_T jBcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 28; j >= 0; j--) {
        jBcol = 12 * j;
        jAcol = 29 * j - 1;
        for (k = j + 2; k < 30; k++) {
            kBcol = (k - 1) * 12;
            if (A[k + jAcol] != 0.0) {
                for (i = 0; i < 12; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[k + jAcol] * B[i + kBcol];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mrdivide_h(const real_T A[348],
                                                                   const real_T B[841],
                                                                   real_T y[348])
{
    int32_T jp;
    int32_T ipiv[29];
    int32_T info;
    real_T temp;
    int32_T xi;
    memcpy(&torqueBalancingYoga_B.b_A[0], &B[0], 841U * sizeof(real_T));
    torqueBalancingYoga_xgetrf_i(torqueBalancingYoga_B.b_A, ipiv, &info);
    memcpy(&y[0], &A[0], 348U * sizeof(real_T));
    torqueBalancingYoga_xtrsm_o(torqueBalancingYoga_B.b_A, y);
    torqueBalancingYoga_xtrsm_oz(torqueBalancingYoga_B.b_A, y);
    for (info = 27; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            jp = ipiv[info] - 1;
            for (xi = 0; xi < 12; xi++) {
                temp = y[12 * info + xi];
                y[xi + 12 * info] = y[12 * jp + xi];
                y[xi + 12 * jp] = temp;
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_mrdivide_hg(real_T A[138],
                                                                    const real_T B[36])
{
    real_T b_A[36];
    int32_T ipiv[6];
    int32_T info;
    real_T b_temp;
    int32_T jBcol;
    int32_T jAcol;
    int32_T kBcol;
    int32_T k;
    int32_T i;
    int32_T tmp;
    memcpy(&b_A[0], &B[0], 36U * sizeof(real_T));
    torqueBalancingYoga_xzgetrf(b_A, ipiv, &info);
    for (info = 0; info < 6; info++) {
        jBcol = 23 * info;
        jAcol = 6 * info;
        for (k = 1; k <= info; k++) {
            kBcol = (k - 1) * 23;
            if (b_A[(k + jAcol) - 1] != 0.0) {
                for (i = 0; i < 23; i++) {
                    tmp = i + jBcol;
                    A[tmp] -= b_A[(k + jAcol) - 1] * A[i + kBcol];
                }
            }
        }

        b_temp = 1.0 / b_A[info + jAcol];
        for (jAcol = 0; jAcol < 23; jAcol++) {
            tmp = jAcol + jBcol;
            A[tmp] *= b_temp;
        }
    }

    for (info = 5; info >= 0; info--) {
        jBcol = 23 * info;
        jAcol = 6 * info - 1;
        for (k = info + 2; k < 7; k++) {
            kBcol = (k - 1) * 23;
            if (b_A[k + jAcol] != 0.0) {
                for (i = 0; i < 23; i++) {
                    A[i + jBcol] -= b_A[k + jAcol] * A[i + kBcol];
                }
            }
        }
    }

    for (info = 4; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            jBcol = ipiv[info] - 1;
            for (jAcol = 0; jAcol < 23; jAcol++) {
                b_temp = A[23 * info + jAcol];
                A[jAcol + 23 * info] = A[23 * jBcol + jAcol];
                A[jAcol + 23 * jBcol] = b_temp;
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_pinvDamped(const real_T A[276],
                                                                   real_T regDamp,
                                                                   real_T pinvDampA[276])
{
    real_T b_A[144];
    int8_T ipiv[12];
    int32_T j;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    int32_T jBcol;
    int32_T kBcol;
    static const int8_T b_b[144] = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    for (j = 0; j < 12; j++) {
        for (jBcol = 0; jBcol < 12; jBcol++) {
            smax = 0.0;
            for (iy = 0; iy < 23; iy++) {
                smax += A[12 * iy + j] * A[12 * iy + jBcol];
            }

            b_A[j + 12 * jBcol] = (real_T) b_b[12 * jBcol + j] * regDamp + smax;
        }

        ipiv[j] = (int8_T)(1 + j);
    }

    for (j = 0; j < 11; j++) {
        jBcol = j * 13;
        iy = 0;
        ix = jBcol;
        smax = std::abs(b_A[jBcol]);
        for (kBcol = 2; kBcol <= 12 - j; kBcol++) {
            ix++;
            s = std::abs(b_A[ix]);
            if (s > smax) {
                iy = kBcol - 1;
                smax = s;
            }
        }

        if (b_A[jBcol + iy] != 0.0) {
            if (iy != 0) {
                iy += j;
                ipiv[j] = (int8_T)(iy + 1);
                ix = j;
                for (kBcol = 0; kBcol < 12; kBcol++) {
                    smax = b_A[ix];
                    b_A[ix] = b_A[iy];
                    b_A[iy] = smax;
                    ix += 12;
                    iy += 12;
                }
            }

            iy = (jBcol - j) + 12;
            for (ix = jBcol + 1; ix < iy; ix++) {
                b_A[ix] /= b_A[jBcol];
            }
        }

        iy = jBcol;
        ix = jBcol + 12;
        for (kBcol = 1; kBcol <= 11 - j; kBcol++) {
            smax = b_A[ix];
            if (b_A[ix] != 0.0) {
                c_ix = jBcol + 1;
                d = (iy - j) + 24;
                for (ijA = 13 + iy; ijA < d; ijA++) {
                    b_A[ijA] += b_A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 12;
            iy += 12;
        }
    }

    for (j = 0; j < 12; j++) {
        for (jBcol = 0; jBcol < 23; jBcol++) {
            pinvDampA[jBcol + 23 * j] = A[12 * jBcol + j];
        }
    }

    for (j = 0; j < 12; j++) {
        jBcol = 23 * j;
        iy = 12 * j;
        for (ix = 1; ix <= j; ix++) {
            kBcol = (ix - 1) * 23;
            if (b_A[(ix + iy) - 1] != 0.0) {
                for (c_ix = 0; c_ix < 23; c_ix++) {
                    d = c_ix + jBcol;
                    pinvDampA[d] -= b_A[(ix + iy) - 1] * pinvDampA[c_ix + kBcol];
                }
            }
        }

        smax = 1.0 / b_A[j + iy];
        for (iy = 0; iy < 23; iy++) {
            d = iy + jBcol;
            pinvDampA[d] *= smax;
        }
    }

    for (j = 11; j >= 0; j--) {
        jBcol = 23 * j;
        iy = 12 * j - 1;
        for (ix = j + 2; ix < 13; ix++) {
            kBcol = (ix - 1) * 23;
            if (b_A[ix + iy] != 0.0) {
                for (c_ix = 0; c_ix < 23; c_ix++) {
                    pinvDampA[c_ix + jBcol] -= b_A[ix + iy] * pinvDampA[c_ix + kBcol];
                }
            }
        }
    }

    for (j = 10; j >= 0; j--) {
        if (j + 1 != ipiv[j]) {
            jBcol = ipiv[j] - 1;
            for (iy = 0; iy < 23; iy++) {
                smax = pinvDampA[23 * j + iy];
                pinvDampA[iy + 23 * j] = pinvDampA[23 * jBcol + iy];
                pinvDampA[iy + 23 * jBcol] = smax;
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_eye_a0(real_T I[529])
{
    int32_T k;
    memset(&I[0], 0, 529U * sizeof(real_T));
    for (k = 0; k < 23; k++) {
        I[k + 23 * k] = 1.0;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_blkdiag(const real_T varargin_1[9],
                                                                const real_T varargin_2[9],
                                                                real_T y[36])
{
    int32_T i;
    int32_T y_tmp;
    int32_T y_tmp_0;
    memset(&y[0], 0, 36U * sizeof(real_T));
    for (i = 0; i < 3; i++) {
        y[6 * i] = varargin_1[3 * i];
        y_tmp_0 = 6 * (3 + i);
        y[3 + y_tmp_0] = varargin_2[3 * i];
        y_tmp = 3 * i + 1;
        y[1 + 6 * i] = varargin_1[y_tmp];
        y[4 + y_tmp_0] = varargin_2[y_tmp];
        y_tmp = 3 * i + 2;
        y[2 + 6 * i] = varargin_1[y_tmp];
        y[5 + y_tmp_0] = varargin_2[y_tmp];
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_blkdiag_a(const real_T varargin_1[114],
                                                                  const real_T varargin_2[114],
                                                                  real_T y[456])
{
    int32_T i;
    int32_T i_0;
    int32_T y_tmp;
    memset(&y[0], 0, 456U * sizeof(real_T));
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (i = 0; i < 19; i++) {
            y_tmp = 19 * i_0 + i;
            y[i + 38 * i_0] = varargin_1[y_tmp];
            y[(i + 38 * (6 + i_0)) + 19] = varargin_2[y_tmp];
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_diag(const real_T v[23], real_T d[529])
{
    int32_T j;
    memset(&d[0], 0, 529U * sizeof(real_T));
    for (j = 0; j < 23; j++) {
        d[j + 23 * j] = v[j];
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xnrm2_mo(int32_T n,
                                                                   const real_T x[529],
                                                                   int32_T ix0)
{
    real_T y;
    real_T scale;
    int32_T kend;
    real_T absxk;
    real_T t;
    int32_T k;
    y = 0.0;
    if (!(n < 1)) {
        if (n == 1) {
            y = std::abs(x[ix0 - 1]);
        }
        else {
            scale = 3.3121686421112381E-170;
            kend = (ix0 + n) - 1;
            for (k = ix0; k <= kend; k++) {
                absxk = std::abs(x[k - 1]);
                if (absxk > scale) {
                    t = scale / absxk;
                    y = y * t * t + 1.0;
                    scale = absxk;
                }
                else {
                    t = absxk / scale;
                    y += t * t;
                }
            }

            y = scale * std::sqrt(y);
        }
    }

    return y;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xnrm2_mot(int32_T n,
                                                                    const real_T x[23],
                                                                    int32_T ix0)
{
    real_T y;
    real_T scale;
    int32_T kend;
    real_T absxk;
    real_T t;
    int32_T k;
    y = 0.0;
    if (!(n < 1)) {
        if (n == 1) {
            y = std::abs(x[ix0 - 1]);
        }
        else {
            scale = 3.3121686421112381E-170;
            kend = (ix0 + n) - 1;
            for (k = ix0; k <= kend; k++) {
                absxk = std::abs(x[k - 1]);
                if (absxk > scale) {
                    t = scale / absxk;
                    y = y * t * t + 1.0;
                    scale = absxk;
                }
                else {
                    t = absxk / scale;
                    y += t * t;
                }
            }

            y = scale * std::sqrt(y);
        }
    }

    return y;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYog_xaxpy_ercdwe(int32_T n,
                                                                    real_T a,
                                                                    const real_T x[23],
                                                                    int32_T ix0,
                                                                    real_T y[529],
                                                                    int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * x[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy_ercdw(int32_T n,
                                                                    real_T a,
                                                                    const real_T x[529],
                                                                    int32_T ix0,
                                                                    real_T y[23],
                                                                    int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * x[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
real_T torqueBalancingYogaModelClass::torqueBalancingYoga_xdotc_lf(int32_T n,
                                                                   const real_T x[529],
                                                                   int32_T ix0,
                                                                   const real_T y[529],
                                                                   int32_T iy0)
{
    real_T d;
    int32_T ix;
    int32_T iy;
    int32_T k;
    d = 0.0;
    if (!(n < 1)) {
        ix = ix0;
        iy = iy0;
        for (k = 1; k <= n; k++) {
            d += x[ix - 1] * y[iy - 1];
            ix++;
            iy++;
        }
    }

    return d;
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xaxpy_ercd(int32_T n,
                                                                   real_T a,
                                                                   int32_T ix0,
                                                                   real_T y[529],
                                                                   int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    int32_T k;
    if (!((n < 1) || (a == 0.0))) {
        ix = ix0 - 1;
        iy = iy0 - 1;
        for (k = 0; k < n; k++) {
            y[iy] += a * y[ix];
            ix++;
            iy++;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xscal_mr(real_T a,
                                                                 real_T x[529],
                                                                 int32_T ix0)
{
    int32_T k;
    for (k = ix0; k <= ix0 + 22; k++) {
        x[k - 1] *= a;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xswap_mc0p(real_T x[529],
                                                                   int32_T ix0,
                                                                   int32_T iy0)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 23; k++) {
        temp = x[ix];
        x[ix] = x[iy];
        x[iy] = temp;
        ix++;
        iy++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xrot_gt(real_T x[529],
                                                                int32_T ix0,
                                                                int32_T iy0,
                                                                real_T c,
                                                                real_T s)
{
    int32_T ix;
    int32_T iy;
    real_T temp;
    int32_T k;
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < 23; k++) {
        temp = c * x[ix] + s * x[iy];
        x[iy] = c * x[iy] - s * x[ix];
        x[ix] = temp;
        iy++;
        ix++;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_svd_b(const real_T A[529],
                                                              real_T U[529],
                                                              real_T s[23],
                                                              real_T V[529])
{
    real_T b_s[23];
    real_T e[23];
    real_T work[23];
    int32_T qq;
    boolean_T apply_transform;
    real_T nrm;
    int32_T qjj;
    int32_T qp1jj;
    int32_T m;
    int32_T qp1q;
    real_T rt;
    real_T ztest;
    real_T smm1;
    real_T emm1;
    real_T sqds;
    real_T shift;
    int32_T k_ii;
    int32_T exitg1;
    boolean_T exitg2;
    memcpy(&torqueBalancingYoga_B.b_A_b[0], &A[0], 529U * sizeof(real_T));
    memset(&b_s[0], 0, 23U * sizeof(real_T));
    memset(&e[0], 0, 23U * sizeof(real_T));
    memset(&work[0], 0, 23U * sizeof(real_T));
    memset(&U[0], 0, 529U * sizeof(real_T));
    memset(&torqueBalancingYoga_B.Vf[0], 0, 529U * sizeof(real_T));
    for (m = 0; m < 22; m++) {
        qq = 23 * m + m;
        apply_transform = false;
        nrm = torqueBalancingYoga_xnrm2_mo(23 - m, torqueBalancingYoga_B.b_A_b, qq + 1);
        if (nrm > 0.0) {
            apply_transform = true;
            if (torqueBalancingYoga_B.b_A_b[qq] < 0.0) {
                b_s[m] = -nrm;
            }
            else {
                b_s[m] = nrm;
            }

            if (std::abs(b_s[m]) >= 1.0020841800044864E-292) {
                nrm = 1.0 / b_s[m];
                qp1q = (qq - m) + 23;
                for (qjj = qq; qjj < qp1q; qjj++) {
                    torqueBalancingYoga_B.b_A_b[qjj] *= nrm;
                }
            }
            else {
                qp1q = (qq - m) + 23;
                for (qjj = qq; qjj < qp1q; qjj++) {
                    torqueBalancingYoga_B.b_A_b[qjj] /= b_s[m];
                }
            }

            torqueBalancingYoga_B.b_A_b[qq]++;
            b_s[m] = -b_s[m];
        }
        else {
            b_s[m] = 0.0;
        }

        for (qp1q = m + 1; qp1q + 1 < 24; qp1q++) {
            qjj = 23 * qp1q + m;
            if (apply_transform) {
                torqueBalancingYoga_xaxpy_ercd(
                    23 - m,
                    -(torqueBalancingYoga_xdotc_lf(23 - m,
                                                   torqueBalancingYoga_B.b_A_b,
                                                   qq + 1,
                                                   torqueBalancingYoga_B.b_A_b,
                                                   qjj + 1)
                      / torqueBalancingYoga_B.b_A_b[m + 23 * m]),
                    qq + 1,
                    torqueBalancingYoga_B.b_A_b,
                    qjj + 1);
            }

            e[qp1q] = torqueBalancingYoga_B.b_A_b[qjj];
        }

        for (qq = m; qq + 1 < 24; qq++) {
            U[qq + 23 * m] = torqueBalancingYoga_B.b_A_b[23 * m + qq];
        }

        if (m + 1 <= 21) {
            nrm = torqueBalancingYoga_xnrm2_mot(22 - m, e, m + 2);
            if (nrm == 0.0) {
                e[m] = 0.0;
            }
            else {
                if (e[m + 1] < 0.0) {
                    e[m] = -nrm;
                }
                else {
                    e[m] = nrm;
                }

                nrm = e[m];
                if (std::abs(e[m]) >= 1.0020841800044864E-292) {
                    nrm = 1.0 / e[m];
                    for (qq = m + 1; qq < 23; qq++) {
                        e[qq] *= nrm;
                    }
                }
                else {
                    for (qq = m + 1; qq < 23; qq++) {
                        e[qq] /= nrm;
                    }
                }

                e[m + 1]++;
                e[m] = -e[m];
                for (qq = m + 1; qq + 1 < 24; qq++) {
                    work[qq] = 0.0;
                }

                for (qq = m + 1; qq + 1 < 24; qq++) {
                    torqueBalancingYoga_xaxpy_ercdw(
                        22 - m, e[qq], torqueBalancingYoga_B.b_A_b, (m + 23 * qq) + 2, work, m + 2);
                }

                for (qq = m + 1; qq + 1 < 24; qq++) {
                    torqueBalancingYog_xaxpy_ercdwe(22 - m,
                                                    -e[qq] / e[m + 1],
                                                    work,
                                                    m + 2,
                                                    torqueBalancingYoga_B.b_A_b,
                                                    (m + 23 * qq) + 2);
                }
            }

            for (qq = m + 1; qq + 1 < 24; qq++) {
                torqueBalancingYoga_B.Vf[qq + 23 * m] = e[qq];
            }
        }
    }

    m = 21;
    b_s[22] = torqueBalancingYoga_B.b_A_b[528];
    e[21] = torqueBalancingYoga_B.b_A_b[527];
    e[22] = 0.0;
    memset(&U[506], 0, 23U * sizeof(real_T));
    U[528] = 1.0;
    for (qp1q = 21; qp1q >= 0; qp1q--) {
        qq = 23 * qp1q + qp1q;
        if (b_s[qp1q] != 0.0) {
            for (qp1jj = qp1q + 1; qp1jj + 1 < 24; qp1jj++) {
                qjj = (23 * qp1jj + qp1q) + 1;
                torqueBalancingYoga_xaxpy_ercd(
                    23 - qp1q,
                    -(torqueBalancingYoga_xdotc_lf(23 - qp1q, U, qq + 1, U, qjj) / U[qq]),
                    qq + 1,
                    U,
                    qjj);
            }

            for (qjj = qp1q; qjj + 1 < 24; qjj++) {
                U[qjj + 23 * qp1q] = -U[23 * qp1q + qjj];
            }

            U[qq]++;
            for (qq = 1; qq <= qp1q; qq++) {
                U[(qq + 23 * qp1q) - 1] = 0.0;
            }
        }
        else {
            memset(&U[qp1q * 23], 0, 23U * sizeof(real_T));
            U[qq] = 1.0;
        }
    }

    for (qq = 22; qq >= 0; qq--) {
        if ((qq + 1 <= 21) && (e[qq] != 0.0)) {
            qp1q = (23 * qq + qq) + 2;
            for (qjj = qq + 1; qjj + 1 < 24; qjj++) {
                qp1jj = (23 * qjj + qq) + 2;
                torqueBalancingYoga_xaxpy_ercd(
                    22 - qq,
                    -(torqueBalancingYoga_xdotc_lf(
                          22 - qq, torqueBalancingYoga_B.Vf, qp1q, torqueBalancingYoga_B.Vf, qp1jj)
                      / torqueBalancingYoga_B.Vf[qp1q - 1]),
                    qp1q,
                    torqueBalancingYoga_B.Vf,
                    qp1jj);
            }
        }

        memset(&torqueBalancingYoga_B.Vf[qq * 23], 0, 23U * sizeof(real_T));
        torqueBalancingYoga_B.Vf[qq + 23 * qq] = 1.0;
    }

    for (qq = 0; qq < 23; qq++) {
        ztest = e[qq];
        if (b_s[qq] != 0.0) {
            rt = std::abs(b_s[qq]);
            nrm = b_s[qq] / rt;
            b_s[qq] = rt;
            if (qq + 1 < 23) {
                ztest = e[qq] / nrm;
            }

            torqueBalancingYoga_xscal_mr(nrm, U, 1 + 23 * qq);
        }

        if ((qq + 1 < 23) && (ztest != 0.0)) {
            rt = std::abs(ztest);
            nrm = rt / ztest;
            ztest = rt;
            b_s[qq + 1] *= nrm;
            torqueBalancingYoga_xscal_mr(nrm, torqueBalancingYoga_B.Vf, 1 + 23 * (qq + 1));
        }

        e[qq] = ztest;
    }

    qq = 0;
    nrm = 0.0;
    for (qp1q = 0; qp1q < 23; qp1q++) {
        ztest = std::abs(b_s[qp1q]);
        rt = std::abs(e[qp1q]);
        if (ztest > rt) {
            rt = ztest;
        }

        if (!(nrm > rt)) {
            nrm = rt;
        }
    }

    while ((m + 2 > 0) && (!(qq >= 75))) {
        qp1jj = m + 1;
        do {
            exitg1 = 0;
            qp1q = qp1jj;
            if (qp1jj == 0) {
                exitg1 = 1;
            }
            else {
                rt = std::abs(e[qp1jj - 1]);
                if ((rt
                     <= (std::abs(b_s[qp1jj - 1]) + std::abs(b_s[qp1jj])) * 2.2204460492503131E-16)
                    || ((rt <= 1.0020841800044864E-292)
                        || ((qq > 20) && (rt <= 2.2204460492503131E-16 * nrm)))) {
                    e[qp1jj - 1] = 0.0;
                    exitg1 = 1;
                }
                else {
                    qp1jj--;
                }
            }
        } while (exitg1 == 0);

        if (m + 1 == qp1jj) {
            qp1jj = 4;
        }
        else {
            qjj = m + 2;
            k_ii = m + 2;
            exitg2 = false;
            while ((!exitg2) && (k_ii >= qp1jj)) {
                qjj = k_ii;
                if (k_ii == qp1jj) {
                    exitg2 = true;
                }
                else {
                    rt = 0.0;
                    if (k_ii < m + 2) {
                        rt = std::abs(e[k_ii - 1]);
                    }

                    if (k_ii > qp1jj + 1) {
                        rt += std::abs(e[k_ii - 2]);
                    }

                    ztest = std::abs(b_s[k_ii - 1]);
                    if ((ztest <= 2.2204460492503131E-16 * rt)
                        || (ztest <= 1.0020841800044864E-292)) {
                        b_s[k_ii - 1] = 0.0;
                        exitg2 = true;
                    }
                    else {
                        k_ii--;
                    }
                }
            }

            if (qjj == qp1jj) {
                qp1jj = 3;
            }
            else if (m + 2 == qjj) {
                qp1jj = 1;
            }
            else {
                qp1jj = 2;
                qp1q = qjj;
            }
        }

        switch (qp1jj) {
            case 1:
                rt = e[m];
                e[m] = 0.0;
                for (qjj = m; qjj + 1 >= qp1q + 1; qjj--) {
                    ztest = b_s[qjj];
                    torqueBalancingYoga_xrotg(&ztest, &rt, &sqds, &smm1);
                    b_s[qjj] = ztest;
                    if (qjj + 1 > qp1q + 1) {
                        rt = e[qjj - 1] * -smm1;
                        e[qjj - 1] *= sqds;
                    }

                    torqueBalancingYoga_xrot_gt(
                        torqueBalancingYoga_B.Vf, 1 + 23 * qjj, 1 + 23 * (m + 1), sqds, smm1);
                }
                break;

            case 2:
                rt = e[qp1q - 1];
                e[qp1q - 1] = 0.0;
                for (qjj = qp1q; qjj < m + 2; qjj++) {
                    ztest = b_s[qjj];
                    torqueBalancingYoga_xrotg(&ztest, &rt, &sqds, &smm1);
                    b_s[qjj] = ztest;
                    rt = -smm1 * e[qjj];
                    e[qjj] *= sqds;
                    torqueBalancingYoga_xrot_gt(U, 1 + 23 * qjj, 1 + 23 * (qp1q - 1), sqds, smm1);
                }
                break;

            case 3:
                ztest = std::abs(b_s[m + 1]);
                rt = std::abs(b_s[m]);
                if (ztest > rt) {
                    rt = ztest;
                }

                ztest = std::abs(e[m]);
                if (rt > ztest) {
                    ztest = rt;
                }

                rt = std::abs(b_s[qp1q]);
                if (ztest > rt) {
                    rt = ztest;
                }

                ztest = std::abs(e[qp1q]);
                if (rt > ztest) {
                    ztest = rt;
                }

                rt = b_s[m + 1] / ztest;
                smm1 = b_s[m] / ztest;
                emm1 = e[m] / ztest;
                sqds = b_s[qp1q] / ztest;
                smm1 = ((smm1 + rt) * (smm1 - rt) + emm1 * emm1) / 2.0;
                emm1 *= rt;
                emm1 *= emm1;
                if ((smm1 != 0.0) || (emm1 != 0.0)) {
                    shift = std::sqrt(smm1 * smm1 + emm1);
                    if (smm1 < 0.0) {
                        shift = -shift;
                    }

                    shift = emm1 / (smm1 + shift);
                }
                else {
                    shift = 0.0;
                }

                rt = (sqds + rt) * (sqds - rt) + shift;
                ztest = e[qp1q] / ztest * sqds;
                for (qjj = qp1q + 1; qjj <= m + 1; qjj++) {
                    torqueBalancingYoga_xrotg(&rt, &ztest, &sqds, &smm1);
                    if (qjj > qp1q + 1) {
                        e[qjj - 2] = rt;
                    }

                    rt = b_s[qjj - 1] * sqds + e[qjj - 1] * smm1;
                    e[qjj - 1] = e[qjj - 1] * sqds - b_s[qjj - 1] * smm1;
                    ztest = smm1 * b_s[qjj];
                    b_s[qjj] *= sqds;
                    torqueBalancingYoga_xrot_gt(
                        torqueBalancingYoga_B.Vf, 1 + 23 * (qjj - 1), 1 + 23 * qjj, sqds, smm1);
                    torqueBalancingYoga_xrotg(&rt, &ztest, &sqds, &smm1);
                    b_s[qjj - 1] = rt;
                    rt = e[qjj - 1] * sqds + smm1 * b_s[qjj];
                    b_s[qjj] = e[qjj - 1] * -smm1 + sqds * b_s[qjj];
                    ztest = smm1 * e[qjj];
                    e[qjj] *= sqds;
                    torqueBalancingYoga_xrot_gt(U, 1 + 23 * (qjj - 1), 1 + 23 * qjj, sqds, smm1);
                }

                e[m] = rt;
                qq++;
                break;

            default:
                if (b_s[qp1q] < 0.0) {
                    b_s[qp1q] = -b_s[qp1q];
                    torqueBalancingYoga_xscal_mr(-1.0, torqueBalancingYoga_B.Vf, 1 + 23 * qp1q);
                }

                qq = qp1q + 1;
                while ((qp1q + 1 < 23) && (b_s[qp1q] < b_s[qq])) {
                    rt = b_s[qp1q];
                    b_s[qp1q] = b_s[qq];
                    b_s[qq] = rt;
                    torqueBalancingYoga_xswap_mc0p(
                        torqueBalancingYoga_B.Vf, 1 + 23 * qp1q, 1 + 23 * (qp1q + 1));
                    torqueBalancingYoga_xswap_mc0p(U, 1 + 23 * qp1q, 1 + 23 * (qp1q + 1));
                    qp1q = qq;
                    qq++;
                }

                qq = 0;
                m--;
                break;
        }
    }

    for (m = 0; m < 23; m++) {
        s[m] = b_s[m];
        memcpy(&V[m * 23], &torqueBalancingYoga_B.Vf[m * 23], 23U * sizeof(real_T));
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xgemm(int32_T k,
                                                              const real_T A[529],
                                                              const real_T B[529],
                                                              real_T C[529])
{
    int32_T br;
    int32_T ar;
    int32_T ia;
    int32_T ic;
    int32_T b;
    int32_T ib;
    int32_T b_ic;
    for (br = 0; br <= 507; br += 23) {
        for (ic = br; ic < br + 23; ic++) {
            C[ic] = 0.0;
        }
    }

    br = -1;
    for (ic = 0; ic <= 507; ic += 23) {
        ar = -1;
        br++;
        b = ((k - 1) * 23 + br) + 1;
        for (ib = br; ib + 1 <= b; ib += 23) {
            if (B[ib] != 0.0) {
                ia = ar;
                for (b_ic = ic; b_ic < ic + 23; b_ic++) {
                    ia++;
                    C[b_ic] += B[ib] * A[ia];
                }
            }

            ar += 23;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_pinv_o(const real_T A[529],
                                                               real_T tol,
                                                               real_T X[529])
{
    int32_T r;
    int32_T vcol;
    real_T s[23];
    int32_T j;
    memset(&X[0], 0, 529U * sizeof(real_T));
    torqueBalancingYoga_svd_b(A, torqueBalancingYoga_B.U, s, torqueBalancingYoga_B.V);
    r = 0;
    vcol = 1;
    while ((vcol < 24) && (s[vcol - 1] > tol)) {
        r++;
        vcol++;
    }

    if (r > 0) {
        vcol = 1;
        for (j = 1; j <= r; j++) {
            torqueBalancingYoga_xscal_mr(1.0 / s[j - 1], torqueBalancingYoga_B.V, vcol);
            vcol += 23;
        }

        torqueBalancingYoga_xgemm(r, torqueBalancingYoga_B.V, torqueBalancingYoga_B.U, X);
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xgetrf_iq(real_T A[529],
                                                                  int32_T ipiv[23],
                                                                  int32_T* info)
{
    int32_T j;
    int32_T c;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T k;
    int32_T jA;
    int32_T b_ix;
    int32_T d;
    int32_T ijA;
    for (j = 0; j < 23; j++) {
        ipiv[j] = 1 + j;
    }

    *info = 0;
    for (j = 0; j < 22; j++) {
        c = j * 24;
        jA = 1;
        ix = c;
        smax = std::abs(A[c]);
        for (k = 2; k <= 23 - j; k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                jA = k;
                smax = s;
            }
        }

        if (A[(c + jA) - 1] != 0.0) {
            if (jA - 1 != 0) {
                ipiv[j] = j + jA;
                torqueBalancingYoga_xswap(A, j + 1, j + jA);
            }

            jA = (c - j) + 23;
            for (ix = c + 1; ix < jA; ix++) {
                A[ix] /= A[c];
            }
        }
        else {
            *info = j + 1;
        }

        jA = c;
        ix = c + 23;
        for (k = 1; k <= 22 - j; k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                b_ix = c + 1;
                d = (jA - j) + 46;
                for (ijA = 24 + jA; ijA < d; ijA++) {
                    A[ijA] += A[b_ix] * -smax;
                    b_ix++;
                }
            }

            ix += 23;
            jA += 23;
        }
    }

    if ((*info == 0) && (!(A[528] != 0.0))) {
        *info = 23;
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm_oz1(const real_T A[529],
                                                                  real_T B[276])
{
    real_T temp;
    int32_T jBcol;
    int32_T jAcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 0; j < 23; j++) {
        jBcol = 12 * j;
        jAcol = 23 * j;
        for (k = 1; k <= j; k++) {
            kBcol = (k - 1) * 12;
            if (A[(k + jAcol) - 1] != 0.0) {
                for (i = 0; i < 12; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[(k + jAcol) - 1] * B[i + kBcol];
                }
            }
        }

        temp = 1.0 / A[j + jAcol];
        for (jAcol = 0; jAcol < 12; jAcol++) {
            tmp = jAcol + jBcol;
            B[tmp] *= temp;
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_xtrsm_oz1v(const real_T A[529],
                                                                   real_T B[276])
{
    int32_T jAcol;
    int32_T jBcol;
    int32_T kBcol;
    int32_T j;
    int32_T k;
    int32_T i;
    int32_T tmp;
    for (j = 22; j >= 0; j--) {
        jBcol = 12 * j;
        jAcol = 23 * j - 1;
        for (k = j + 2; k < 24; k++) {
            kBcol = (k - 1) * 12;
            if (A[k + jAcol] != 0.0) {
                for (i = 0; i < 12; i++) {
                    tmp = i + jBcol;
                    B[tmp] -= A[k + jAcol] * B[i + kBcol];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S6>/Balancing Controller ' */
void torqueBalancingYogaModelClass::torqueBalancingYog_pinvDamped_a(const real_T A[276],
                                                                    real_T regDamp,
                                                                    real_T pinvDampA[276])
{
    int32_T jp;
    int32_T ipiv[23];
    int32_T info;
    real_T temp;
    int32_T xi;
    static const int8_T b[529] = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    for (info = 0; info < 23; info++) {
        for (jp = 0; jp < 23; jp++) {
            temp = 0.0;
            for (xi = 0; xi < 12; xi++) {
                temp += A[23 * xi + info] * A[23 * xi + jp];
            }

            torqueBalancingYoga_B.b_A_c[info + 23 * jp] =
                (real_T) b[23 * jp + info] * regDamp + temp;
        }
    }

    torqueBalancingYoga_xgetrf_iq(torqueBalancingYoga_B.b_A_c, ipiv, &info);
    for (info = 0; info < 23; info++) {
        for (jp = 0; jp < 12; jp++) {
            pinvDampA[jp + 12 * info] = A[23 * jp + info];
        }
    }

    torqueBalancingYoga_xtrsm_oz1(torqueBalancingYoga_B.b_A_c, pinvDampA);
    torqueBalancingYoga_xtrsm_oz1v(torqueBalancingYoga_B.b_A_c, pinvDampA);
    for (info = 21; info >= 0; info--) {
        if (info + 1 != ipiv[info]) {
            jp = ipiv[info] - 1;
            for (xi = 0; xi < 12; xi++) {
                temp = pinvDampA[12 * info + xi];
                pinvDampA[xi + 12 * info] = pinvDampA[12 * jp + xi];
                pinvDampA[xi + 12 * jp] = temp;
            }
        }
    }
}

/* Function for MATLAB Function: '<S144>/Analytical Solution One Foot (unconstrained)' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_invNxN(const real_T x[36], real_T y[36])
{
    int8_T p[6];
    real_T A[36];
    int8_T ipiv[6];
    int32_T b_j;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T b_k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    int32_T pipk;
    for (b_k = 0; b_k < 36; b_k++) {
        y[b_k] = 0.0;
        A[b_k] = x[b_k];
    }

    for (b_k = 0; b_k < 6; b_k++) {
        ipiv[b_k] = (int8_T)(1 + b_k);
    }

    for (b_j = 0; b_j < 5; b_j++) {
        pipk = b_j * 7;
        iy = 0;
        ix = pipk;
        smax = std::abs(A[pipk]);
        for (b_k = 2; b_k <= 6 - b_j; b_k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                iy = b_k - 1;
                smax = s;
            }
        }

        if (A[pipk + iy] != 0.0) {
            if (iy != 0) {
                iy += b_j;
                ipiv[b_j] = (int8_T)(iy + 1);
                ix = b_j;
                for (b_k = 0; b_k < 6; b_k++) {
                    smax = A[ix];
                    A[ix] = A[iy];
                    A[iy] = smax;
                    ix += 6;
                    iy += 6;
                }
            }

            iy = (pipk - b_j) + 6;
            for (ix = pipk + 1; ix < iy; ix++) {
                A[ix] /= A[pipk];
            }
        }

        iy = pipk;
        ix = pipk + 6;
        for (b_k = 1; b_k <= 5 - b_j; b_k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                c_ix = pipk + 1;
                d = (iy - b_j) + 12;
                for (ijA = 7 + iy; ijA < d; ijA++) {
                    A[ijA] += A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 6;
            iy += 6;
        }
    }

    for (b_k = 0; b_k < 6; b_k++) {
        p[b_k] = (int8_T)(1 + b_k);
    }

    for (b_j = 0; b_j < 5; b_j++) {
        if (ipiv[b_j] > 1 + b_j) {
            b_k = ipiv[b_j] - 1;
            pipk = p[b_k];
            p[b_k] = p[b_j];
            p[b_j] = (int8_T) pipk;
        }
    }

    for (b_j = 0; b_j < 6; b_j++) {
        b_k = p[b_j] - 1;
        y[b_j + 6 * b_k] = 1.0;
        for (iy = b_j; iy + 1 < 7; iy++) {
            if (y[6 * b_k + iy] != 0.0) {
                for (ix = iy + 1; ix + 1 < 7; ix++) {
                    y[ix + 6 * b_k] -= y[6 * b_k + iy] * A[6 * iy + ix];
                }
            }
        }
    }

    for (b_j = 0; b_j < 6; b_j++) {
        pipk = 6 * b_j;
        for (iy = 5; iy >= 0; iy--) {
            ix = 6 * iy;
            b_k = iy + pipk;
            if (y[b_k] != 0.0) {
                y[b_k] = y[iy + pipk] / A[iy + ix];
                for (b_k = 0; b_k < iy; b_k++) {
                    c_ix = b_k + pipk;
                    y[c_ix] -= y[iy + pipk] * A[b_k + ix];
                }
            }
        }
    }
}

/* Function for MATLAB Function: '<S146>/Analytical Solution Two Feet (unconstrained)' */
void torqueBalancingYogaModelClass::torqueBalancingYoga_invNxN_l(const real_T x[144], real_T y[144])
{
    int8_T p[12];
    real_T A[144];
    int8_T ipiv[12];
    int32_T b_j;
    int32_T ix;
    real_T smax;
    real_T s;
    int32_T b_k;
    int32_T iy;
    int32_T c_ix;
    int32_T d;
    int32_T ijA;
    int32_T pipk;
    for (b_k = 0; b_k < 144; b_k++) {
        y[b_k] = 0.0;
        A[b_k] = x[b_k];
    }

    for (b_k = 0; b_k < 12; b_k++) {
        ipiv[b_k] = (int8_T)(1 + b_k);
    }

    for (b_j = 0; b_j < 11; b_j++) {
        pipk = b_j * 13;
        iy = 0;
        ix = pipk;
        smax = std::abs(A[pipk]);
        for (b_k = 2; b_k <= 12 - b_j; b_k++) {
            ix++;
            s = std::abs(A[ix]);
            if (s > smax) {
                iy = b_k - 1;
                smax = s;
            }
        }

        if (A[pipk + iy] != 0.0) {
            if (iy != 0) {
                iy += b_j;
                ipiv[b_j] = (int8_T)(iy + 1);
                ix = b_j;
                for (b_k = 0; b_k < 12; b_k++) {
                    smax = A[ix];
                    A[ix] = A[iy];
                    A[iy] = smax;
                    ix += 12;
                    iy += 12;
                }
            }

            iy = (pipk - b_j) + 12;
            for (ix = pipk + 1; ix < iy; ix++) {
                A[ix] /= A[pipk];
            }
        }

        iy = pipk;
        ix = pipk + 12;
        for (b_k = 1; b_k <= 11 - b_j; b_k++) {
            smax = A[ix];
            if (A[ix] != 0.0) {
                c_ix = pipk + 1;
                d = (iy - b_j) + 24;
                for (ijA = 13 + iy; ijA < d; ijA++) {
                    A[ijA] += A[c_ix] * -smax;
                    c_ix++;
                }
            }

            ix += 12;
            iy += 12;
        }
    }

    for (b_k = 0; b_k < 12; b_k++) {
        p[b_k] = (int8_T)(1 + b_k);
    }

    for (b_j = 0; b_j < 11; b_j++) {
        if (ipiv[b_j] > 1 + b_j) {
            b_k = ipiv[b_j] - 1;
            pipk = p[b_k];
            p[b_k] = p[b_j];
            p[b_j] = (int8_T) pipk;
        }
    }

    for (b_j = 0; b_j < 12; b_j++) {
        b_k = p[b_j] - 1;
        y[b_j + 12 * b_k] = 1.0;
        for (iy = b_j; iy + 1 < 13; iy++) {
            if (y[12 * b_k + iy] != 0.0) {
                for (ix = iy + 1; ix + 1 < 13; ix++) {
                    y[ix + 12 * b_k] -= y[12 * b_k + iy] * A[12 * iy + ix];
                }
            }
        }
    }

    for (b_j = 0; b_j < 12; b_j++) {
        pipk = 12 * b_j;
        for (iy = 11; iy >= 0; iy--) {
            ix = 12 * iy;
            b_k = iy + pipk;
            if (y[b_k] != 0.0) {
                y[b_k] = y[iy + pipk] / A[iy + ix];
                for (b_k = 0; b_k < iy; b_k++) {
                    c_ix = b_k + pipk;
                    y[c_ix] -= y[iy + pipk] * A[b_k + ix];
                }
            }
        }
    }
}

/* Model step function */
void torqueBalancingYogaModelClass::step()
{
    /* local block i/o variables */
    real_T rtb_imu_H_link[16];
    real_T rtb_link_H_root[16];
    real_T rtb_Switch[3];
    real_T rtb_imu_H_link_f[16];
    real_T rtb_link_H_root_b[16];
    real_T rtb_Switch_l[3];
    real_T rtb_Switch_c[23];
    real_T rtb_imu_H_link_k[16];
    real_T rtb_link_H_root_n[16];
    real_T rtb_Switch_o[3];
    real_T rtb_CoM6DCoMXYZ2_i[3];
    real_T rtb_imu_H_link_l[16];
    real_T rtb_link_H_root_o[16];
    real_T rtb_Switch_m[3];
    real_T rtb_imu_H_link_a[16];
    real_T rtb_link_H_root_l[16];
    real_T rtb_Switch_j[3];
    real_T rtb_CoM6DCoMXYZ2_b[3];
    boolean_T rtb_Compare;
    boolean_T rtb_Compare_c;
    boolean_T rtb_Compare_i;
    boolean_T rtb_Compare_h;
    boolean_T rtb_Compare_a;
    boolean_T rtb_Compare_l;
    boolean_T rtb_Compare_b;
    boolean_T rtb_Compare_j;
    boolean_T rtb_Compare_f;
    boolean_T rtb_Compare_d;
    boolean_T rtb_Compare_bw;
    boolean_T rtb_Compare_g;
    boolean_T rtb_Compare_p;
    boolean_T rtb_Compare_f4;
    boolean_T rtb_Compare_fn;
    static const real_T b[529] = {
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

    real_T pinvJb[36];
    static const int8_T b_0[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};

    real_T gravityWrench[6];
    real_T AR[36];
    real_T pinvA[72];
    real_T JcMinv[348];
    real_T JcMinvSt[276];
    real_T Pinv_JcMinvSt[276];
    real_T constraintMatrixRightFoot[114];
    real_T ConstraintsMatrix2Feet[456];
    real_T LDotDes[6];
    real_T a[72];
    real_T c_a[138];
    static const int8_T c[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    static const int8_T d[18] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    static const int8_T c_b[667] = {
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    static const int8_T f_a[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};

    static const int8_T g_a[144] = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    static const int8_T d_b[529] = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    real_T fixed_link_CoMDes[4];
    real_T Jc[348];
    real_T pinvJb_0[72];
    static const int8_T b_1[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};

    static const int8_T b_2[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};

    boolean_T res[23];
    real_T rtb_fLdotDesC1C2[12];
    real_T rtb_Switch3[9];
    real_T rtb_Clock1;
    real_T rtb_Switch_jf[16];
    real_T rtb_impedances[23];
    real_T rtb_NA[144];
    real_T rtb_Sigma[276];
    int32_T i;
    real_T tmp[9];
    real_T tmp_0[9];
    real_T tmp_1[2];
    int32_T i_0;
    real_T pinvJb_1[72];
    real_T pinvJb_2[138];
    real_T pinvJb_3[36];
    real_T tmp_2[72];
    real_T tmp_3[144];
    real_T JcMinv_0[12];
    real_T tmp_4[23];
    real_T Pinv_JcMinvSt_0[23];
    real_T invTGamma[23];
    real_T JcMinv_1[144];
    real_T Pinv_JcMinvSt_1[276];
    real_T LDotDes_0[6];
    real_T pinvA_0[12];
    real_T rtb_impedances_0[23];
    real_T tmp_5[38];
    real_T ConstraintsMatrix2Feet_0[38];
    real_T rtb_KDCoM_idx_2;
    real_T rtb_KDCoM_idx_1;
    real_T rtb_KDCoM_idx_0;
    real_T rtb_KPCoM_idx_2;
    real_T rtb_KPCoM_idx_1;
    real_T rtb_KPCoM_idx_0;
    real_T acc_CoM_des_idx_0;
    real_T acc_CoM_des_idx_1;
    real_T acc_CoM_des_idx_2;
    int32_T Gamma_tmp;
    int32_T AR_tmp;
    real_T unusedExpr[276];

    /* MATLAB Function: '<S111>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}' */
    /* MATLAB Function 'controller_QP/Compute joint torques (motor reflected
     * inertia)/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}': '<S154>:1' */
    /* '<S154>:1:3' */
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            i = AR_tmp + 23 * i_0;
            torqueBalancingYoga_B.Gamma_m[i] = 0.0;
            for (Gamma_tmp = 0; Gamma_tmp < 23; Gamma_tmp++) {
                torqueBalancingYoga_B.Gamma_m[i] =
                    torqueBalancingYoga_P.Config.T[23 * Gamma_tmp + AR_tmp]
                        * torqueBalancingYoga_P.Config.Gamma[23 * i_0 + Gamma_tmp]
                    + torqueBalancingYoga_B.Gamma_m[23 * i_0 + AR_tmp];
            }
        }
    }

    torqueBalancingYoga_mrdivide(b, torqueBalancingYoga_B.Gamma_m, torqueBalancingYoga_B.invTGamma);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            i = i_0 + 23 * AR_tmp;
            torqueBalancingYoga_B.Gamma_m[i] = 0.0;
            for (Gamma_tmp = 0; Gamma_tmp < 23; Gamma_tmp++) {
                torqueBalancingYoga_B.Gamma_m[i] =
                    torqueBalancingYoga_P.Config.T[23 * Gamma_tmp + AR_tmp]
                        * torqueBalancingYoga_P.Config.Gamma[23 * i_0 + Gamma_tmp]
                    + torqueBalancingYoga_B.Gamma_m[23 * AR_tmp + i_0];
            }
        }
    }

    torqueBalancingYoga_mrdivide(b, torqueBalancingYoga_B.Gamma_m, torqueBalancingYoga_B.Gamma);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            Gamma_tmp = i_0 + 23 * AR_tmp;
            torqueBalancingYoga_B.Gamma_m[Gamma_tmp] = 0.0;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.Gamma_m[Gamma_tmp] =
                    torqueBalancingYoga_B.Gamma[23 * i + i_0]
                        * torqueBalancingYoga_P.Config.I_m[23 * AR_tmp + i]
                    + torqueBalancingYoga_B.Gamma_m[23 * AR_tmp + i_0];
            }
        }

        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            Gamma_tmp = i_0 + 23 * AR_tmp;
            torqueBalancingYoga_B.reflectedInertia[Gamma_tmp] = 0.0;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.reflectedInertia[Gamma_tmp] =
                    torqueBalancingYoga_B.Gamma_m[23 * i + i_0]
                        * torqueBalancingYoga_B.invTGamma[23 * AR_tmp + i]
                    + torqueBalancingYoga_B.reflectedInertia[23 * AR_tmp + i_0];
            }
        }
    }

    /* End of MATLAB Function: '<S111>/(transpose(T*Gamma))^{-1}*I_m*(T*Gamma)^{-1}' */

    /* S-Function (WBToolbox): '<S156>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S156>/S-Function

    /* S-Function (WBToolbox): '<S31>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S31>/S-Function

    /* RelationalOperator: '<S30>/Compare' incorporates:
     *  Constant: '<S30>/Constant'
     *  Constant: '<S4>/Constant'
     *  S-Function (sfix_bitop): '<S4>/Coordinator'
     */
    torqueBalancingYoga_B.Compare =
        ((torqueBalancingYoga_P.Constant_Value_it & torqueBalancingYoga_P.Coordinator_BitMask)
         != torqueBalancingYoga_P.Constant_Value_n);

    /* RelationalOperator: '<S29>/Compare' incorporates:
     *  Constant: '<S29>/Constant'
     *  Constant: '<S4>/Constant'
     *  S-Function (sfix_bitop): '<S4>/Yoga'
     */
    torqueBalancingYoga_B.Compare_c =
        ((torqueBalancingYoga_P.Constant_Value_it & torqueBalancingYoga_P.Yoga_BitMask)
         != torqueBalancingYoga_P.Constant_Value_dx);

    /* Outputs for Enabled SubSystem: '<S4>/State Machine Yoga' incorporates:
     *  EnablePort: '<S34>/Enable'
     */
    if (torqueBalancingYoga_B.Compare_c) {
        if (!torqueBalancingYoga_DW.StateMachineYoga_MODE) {
            torqueBalancingYoga_DW.StateMachineYoga_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.StateMachineYoga_MODE) {
            torqueBalancingYoga_DW.StateMachineYoga_MODE = false;
        }
    }

    if (torqueBalancingYoga_DW.StateMachineYoga_MODE) {
        /* Clock: '<S34>/Clock1' */
        rtb_Clock1 = (&torqueBalancingYoga_M)->Timing.t[0];

        /* S-Function (WBToolbox): '<S34>/right_foot_wrench' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S34>/right_foot_wrench

        /* S-Function (WBToolbox): '<S34>/left_foot_wrench' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S34>/left_foot_wrench

        /* Inport: '<S34>/jointAngles' */
        memcpy(&torqueBalancingYoga_B.jointAngles[0],
               &torqueBalancingYoga_B.SFunction_d[0],
               23U * sizeof(real_T));

        /* S-Function (WBToolbox): '<S82>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S82>/S-Function

        /* S-Function (WBToolbox): '<S79>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S79>/S-Function

        /* Product: '<S78>/Inv\+' */
        rt_mldivided4x4(torqueBalancingYoga_B.SFunction_k,
                        torqueBalancingYoga_B.SFunction_ao,
                        rtb_imu_H_link_k);

        /* RelationalOperator: '<S88>/Compare' incorporates:
         *  Clock: '<S85>/Clock'
         *  Constant: '<S88>/Constant'
         */
        rtb_Compare_l = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_p);

        /* MATLAB Function: '<S85>/MATLAB Function' */
        torqueBalancingY_MATLABFunction(rtb_imu_H_link_k,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_i,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_i);

        /* S-Function (WBToolbox): '<S83>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S83>/S-Function

        /* Product: '<S78>/Inv\*  ' */
        rt_mldivided4x4(torqueBalancingYoga_B.SFunction_ao,
                        torqueBalancingYoga_B.SFunction_nj,
                        rtb_link_H_root_n);

        /* S-Function (WBToolbox): '<S34>/inertial' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S34>/inertial

        /* RelationalOperator: '<S90>/Compare' incorporates:
         *  Clock: '<S86>/Clock'
         *  Constant: '<S90>/Constant'
         */
        rtb_Compare_b = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_od);

        /* MATLAB Function: '<S86>/MATLAB Function' */
        torqueBalancin_MATLABFunction_l(torqueBalancingYoga_B.inertial_n,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_oa,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_oa);

        /* S-Function (WBToolbox): '<S78>/Neck Position' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S78>/Neck Position

        /* Switch: '<S87>/Switch' incorporates:
         *  Constant: '<S87>/Constant'
         *  Constant: '<S87>/USE_IMU4EST_BASE1'
         *  Gain: '<S87>/Gain'
         */
        if (torqueBalancingYoga_P.Config.CORRECT_NECK_IMU) {
            rtb_Switch_o[0] =
                torqueBalancingYoga_P.Gain_Gain_a * torqueBalancingYoga_B.NeckPosition_p[0];
            rtb_Switch_o[1] =
                torqueBalancingYoga_P.Gain_Gain_a * torqueBalancingYoga_B.NeckPosition_p[1];
            rtb_Switch_o[2] =
                torqueBalancingYoga_P.Gain_Gain_a * torqueBalancingYoga_B.NeckPosition_p[2];
        }
        else {
            rtb_Switch_o[0] = torqueBalancingYoga_P.Constant_Value_c[0];
            rtb_Switch_o[1] = torqueBalancingYoga_P.Constant_Value_c[1];
            rtb_Switch_o[2] = torqueBalancingYoga_P.Constant_Value_c[2];
        }

        /* End of Switch: '<S87>/Switch' */

        /* MATLAB Function: '<S78>/fromImuToHomogeousTransformFCN' */
        fromImuToHomogeousTransformFCN(rtb_imu_H_link_k,
                                       torqueBalancingYoga_B.sf_MATLABFunction_i.s0,
                                       rtb_link_H_root_n,
                                       torqueBalancingYoga_B.sf_MATLABFunction_oa.s0,
                                       torqueBalancingYoga_B.inertial_n,
                                       rtb_Switch_o,
                                       &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfo_fb,
                                       &torqueBalancingYoga_P);

        /* Switch: '<S78>/Switch6' incorporates:
         *  Constant: '<S78>/USE_IMU4EST_BASE'
         */
        if (torqueBalancingYoga_P.Config.USE_IMU4EST_BASE) {
            memcpy(&torqueBalancingYoga_B.Switch6[0],
                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfo_fb.w_H_b[0],
                   sizeof(real_T) << 4U);
        }
        else {
            memcpy(&torqueBalancingYoga_B.Switch6[0], &rtb_link_H_root_n[0], sizeof(real_T) << 4U);
        }

        /* End of Switch: '<S78>/Switch6' */

        /* S-Function (WBToolbox): '<S106>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S106>/S-Function

        /* Selector: '<S71>/CoM6D ->  CoMXYZ2' */
        rtb_CoM6DCoMXYZ2_i[0] = torqueBalancingYoga_B.SFunction_jj[12];
        rtb_CoM6DCoMXYZ2_i[1] = torqueBalancingYoga_B.SFunction_jj[13];
        rtb_CoM6DCoMXYZ2_i[2] = torqueBalancingYoga_B.SFunction_jj[14];

        /* RelationalOperator: '<S102>/Compare' incorporates:
         *  Clock: '<S68>/Clock'
         *  Constant: '<S102>/Constant'
         */
        rtb_Compare_j = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_ph);

        /* MATLAB Function: '<S68>/MATLAB Function' */
        torqueBalancin_MATLABFunction_a(rtb_CoM6DCoMXYZ2_i,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_f,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_f);

        /* RelationalOperator: '<S104>/Compare' incorporates:
         *  Clock: '<S69>/Clock'
         *  Constant: '<S104>/Constant'
         */
        rtb_Compare_f = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_pf);

        /* MATLAB Function: '<S69>/MATLAB Function' */
        torqueBalancin_MATLABFunction_k(torqueBalancingYoga_B.jointAngles,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_d,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_d);

        /* S-Function (WBToolbox): '<S92>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S92>/S-Function

        /* S-Function (WBToolbox): '<S81>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S81>/S-Function

        /* Product: '<S80>/Inv\+' */
        rt_mldivided4x4(torqueBalancingYoga_B.SFunction_oz,
                        torqueBalancingYoga_B.SFunction_i,
                        rtb_imu_H_link_l);

        /* RelationalOperator: '<S98>/Compare' incorporates:
         *  Clock: '<S95>/Clock'
         *  Constant: '<S98>/Constant'
         */
        rtb_Compare_d = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_f);

        /* MATLAB Function: '<S95>/MATLAB Function' */
        torqueBalancingY_MATLABFunction(rtb_imu_H_link_l,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_ad,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_ad);

        /* S-Function (WBToolbox): '<S93>/S-Function' incorporates:
         *  Constant: '<S66>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S93>/S-Function

        /* Product: '<S80>/Inv\*  ' */
        rt_mldivided4x4(torqueBalancingYoga_B.SFunction_i,
                        torqueBalancingYoga_B.SFunction_ef,
                        rtb_link_H_root_o);

        /* RelationalOperator: '<S100>/Compare' incorporates:
         *  Clock: '<S96>/Clock'
         *  Constant: '<S100>/Constant'
         */
        rtb_Compare_bw = ((&torqueBalancingYoga_M)->Timing.t[0]
                          == torqueBalancingYoga_P.CompareToConstant_const_e);

        /* MATLAB Function: '<S96>/MATLAB Function' */
        torqueBalancin_MATLABFunction_l(torqueBalancingYoga_B.inertial_n,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_j,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_j);

        /* S-Function (WBToolbox): '<S80>/Neck Position' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S80>/Neck Position

        /* Switch: '<S97>/Switch' incorporates:
         *  Constant: '<S97>/Constant'
         *  Constant: '<S97>/USE_IMU4EST_BASE1'
         *  Gain: '<S97>/Gain'
         */
        if (torqueBalancingYoga_P.Config.CORRECT_NECK_IMU) {
            rtb_Switch_m[0] =
                torqueBalancingYoga_P.Gain_Gain_h * torqueBalancingYoga_B.NeckPosition_k[0];
            rtb_Switch_m[1] =
                torqueBalancingYoga_P.Gain_Gain_h * torqueBalancingYoga_B.NeckPosition_k[1];
            rtb_Switch_m[2] =
                torqueBalancingYoga_P.Gain_Gain_h * torqueBalancingYoga_B.NeckPosition_k[2];
        }
        else {
            rtb_Switch_m[0] = torqueBalancingYoga_P.Constant_Value_cj[0];
            rtb_Switch_m[1] = torqueBalancingYoga_P.Constant_Value_cj[1];
            rtb_Switch_m[2] = torqueBalancingYoga_P.Constant_Value_cj[2];
        }

        /* End of Switch: '<S97>/Switch' */

        /* MATLAB Function: '<S80>/fromImuToHomogeousTransformFCN' */
        fromImuToHomogeousTransformFCN(rtb_imu_H_link_l,
                                       torqueBalancingYoga_B.sf_MATLABFunction_ad.s0,
                                       rtb_link_H_root_o,
                                       torqueBalancingYoga_B.sf_MATLABFunction_j.s0,
                                       torqueBalancingYoga_B.inertial_n,
                                       rtb_Switch_m,
                                       &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_g,
                                       &torqueBalancingYoga_P);

        /* Switch: '<S80>/Switch6' incorporates:
         *  Constant: '<S80>/USE_IMU4EST_BASE'
         */
        if (torqueBalancingYoga_P.Config.USE_IMU4EST_BASE) {
            memcpy(&torqueBalancingYoga_B.Switch6_g[0],
                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_g.w_H_b[0],
                   sizeof(real_T) << 4U);
        }
        else {
            memcpy(
                &torqueBalancingYoga_B.Switch6_g[0], &rtb_link_H_root_o[0], sizeof(real_T) << 4U);
        }

        /* End of Switch: '<S80>/Switch6' */

        /* S-Function (WBToolbox): '<S107>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S107>/S-Function

        /* MATLAB Function: '<S34>/stateMachineYogaFCN' incorporates:
         *  Selector: '<S72>/CoM6D ->  CoMXYZ2'
         */
        /* MATLAB Function 'Robot State and References/State Machine Yoga/stateMachineYogaFCN':
         * '<S70>:1' */
        /* '<S70>:1:4' */
        if ((!torqueBalancingYoga_DW.state_not_empty) || (!torqueBalancingYoga_DW.tSwitch_not_empty)
            || (!torqueBalancingYoga_DW.w_H_fixedLink_not_empty)
            || (!torqueBalancingYoga_DW.secondYoga_not_empty)) {
            torqueBalancingYoga_DW.state = torqueBalancingYoga_P.Sm.stateAt0;
            torqueBalancingYoga_DW.state_not_empty = true;
            torqueBalancingYoga_DW.tSwitch = 0.0;
            torqueBalancingYoga_DW.tSwitch_not_empty = true;
            torqueBalancingYoga_eye(torqueBalancingYoga_DW.w_H_fixedLink);
            torqueBalancingYoga_DW.w_H_fixedLink_not_empty = true;
            torqueBalancingYoga_DW.secondYoga = false;
            torqueBalancingYoga_DW.secondYoga_not_empty = true;
        }

        torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_B.sf_MATLABFunction_f.s0[0];
        torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_B.sf_MATLABFunction_f.s0[1];
        torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
        torqueBalancingYoga_B.constraints[0] = 1.0;
        torqueBalancingYoga_B.constraints[1] = 1.0;
        torqueBalancingYoga_eye(torqueBalancingYoga_B.w_H_b);
        memcpy(&torqueBalancingYoga_B.qj_des[0],
               &torqueBalancingYoga_B.sf_MATLABFunction_d.s0[0],
               23U * sizeof(real_T));
        for (i = 0; i < 23; i++) {
            rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[13 * i];
        }

        rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[0];
        rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[0];
        rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[13];
        rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[13];
        rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[26];
        rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[26];
        if (torqueBalancingYoga_DW.state == 1.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.tBalancing)
                && (!torqueBalancingYoga_P.Sm.demoOnlyBalancing)) {
                torqueBalancingYoga_DW.state = 2.0;
                if (torqueBalancingYoga_P.Sm.demoStartsOnRightSupport) {
                    for (i_0 = 0; i_0 < 4; i_0++) {
                        for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] = 0.0;
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[i_0 << 2]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                        }
                    }

                    for (i_0 = 0; i_0 < 4; i_0++) {
                        torqueBalancingYoga_DW.w_H_fixedLink[i_0 << 2] = rtb_Switch_jf[i_0 << 2];
                        torqueBalancingYoga_DW.w_H_fixedLink[1 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 1];
                        torqueBalancingYoga_DW.w_H_fixedLink[2 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 2];
                        torqueBalancingYoga_DW.w_H_fixedLink[3 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 3];
                    }

                    torqueBalancingYoga_mrdivide_o(torqueBalancingYoga_DW.w_H_fixedLink,
                                                   torqueBalancingYoga_B.Switch6_g);
                    torqueBalancingYoga_DW.state = 8.0;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 2.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                rtb_impedances[AR_tmp] =
                    torqueBalancingYoga_P.Gain.impedances[(13 * AR_tmp + i_0) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            fixed_link_CoMDes[0] = torqueBalancingYoga_B.CoM_des[0];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            fixed_link_CoMDes[1] = torqueBalancingYoga_B.CoM_des[1];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            fixed_link_CoMDes[2] = torqueBalancingYoga_B.CoM_des[2];
            fixed_link_CoMDes[3] = 1.0;
            torqueBalancingYoga_mldivide_e(torqueBalancingYoga_DW.w_H_fixedLink, fixed_link_CoMDes);
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                torqueBalancingYoga_B.qj_des[AR_tmp] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * AR_tmp + i_0) - 1];
            }

            if ((std::abs(fixed_link_CoMDes[1] - rtb_CoM6DCoMXYZ2_i[1])
                 < torqueBalancingYoga_P.Sm.CoM_threshold)
                && (torqueBalancingYoga_B.wR_WBDT[2]
                    < torqueBalancingYoga_P.Sm.wrench_thresholdContactOff)) {
                torqueBalancingYoga_DW.state = 3.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
            }
        }

        if (torqueBalancingYoga_DW.state == 3.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 0.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            if (rtb_Clock1
                > torqueBalancingYoga_DW.tSwitch + torqueBalancingYoga_P.Sm.tBalancingBeforeYoga) {
                torqueBalancingYoga_DW.state = 4.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                if (torqueBalancingYoga_P.Sm.skipYoga) {
                    torqueBalancingYoga_DW.state = 5.0;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 4.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 0.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            for (i = 0; i < 25; i++) {
                if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_leftYogaRef[i]
                                      + torqueBalancingYoga_DW.tSwitch)
                    && (rtb_Clock1 <= torqueBalancingYoga_P.Sm.joints_leftYogaRef[i + 1]
                                          + torqueBalancingYoga_DW.tSwitch)
                    && (!torqueBalancingYoga_DW.secondYoga)) {
                    for (i_0 = 0; i_0 < 23; i_0++) {
                        torqueBalancingYoga_B.qj_des[i_0] =
                            torqueBalancingYoga_P.Sm.joints_leftYogaRef[(1 + i_0) * 26 + i];
                    }
                }
                else {
                    if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_leftSecondYogaRef[i]
                                          + torqueBalancingYoga_DW.tSwitch)
                        && (rtb_Clock1 <= torqueBalancingYoga_P.Sm.joints_leftSecondYogaRef[i + 1]
                                              + torqueBalancingYoga_DW.tSwitch)
                        && torqueBalancingYoga_DW.secondYoga) {
                        for (i_0 = 0; i_0 < 23; i_0++) {
                            torqueBalancingYoga_B.qj_des[i_0] =
                                torqueBalancingYoga_P.Sm.joints_leftYogaRef[(1 + i_0) * 26 + i];
                        }
                    }
                }
            }

            if ((rtb_Clock1
                 > torqueBalancingYoga_P.Sm.joints_leftYogaRef[25] + torqueBalancingYoga_DW.tSwitch)
                && (!torqueBalancingYoga_DW.secondYoga)) {
                for (i_0 = 0; i_0 < 23; i_0++) {
                    torqueBalancingYoga_B.qj_des[i_0] =
                        torqueBalancingYoga_P.Sm.joints_leftYogaRef[(1 + i_0) * 26 + 25];
                }

                if ((rtb_Clock1
                     > ((torqueBalancingYoga_P.Sm.joints_leftYogaRef[25]
                         + torqueBalancingYoga_DW.tSwitch)
                        + torqueBalancingYoga_P.Sm
                              .smoothingTimeCoM_Joints[(int32_T) torqueBalancingYoga_DW.state - 1])
                           + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves)
                    && torqueBalancingYoga_P.Sm.repeatYogaMoveset) {
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                    torqueBalancingYoga_DW.secondYoga = true;
                }
                else {
                    if ((rtb_Clock1 > ((torqueBalancingYoga_P.Sm.joints_leftYogaRef[25]
                                        + torqueBalancingYoga_DW.tSwitch)
                                       + torqueBalancingYoga_P.Sm.smoothingTimeCoM_Joints
                                             [(int32_T) torqueBalancingYoga_DW.state - 1])
                                          + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves)
                        && (!torqueBalancingYoga_P.Sm.repeatYogaMoveset)) {
                        torqueBalancingYoga_DW.state = 5.0;
                        torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                        torqueBalancingYoga_DW.secondYoga = false;
                    }
                }
            }

            if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_leftSecondYogaRef[25]
                                  + torqueBalancingYoga_DW.tSwitch)
                && torqueBalancingYoga_DW.secondYoga) {
                for (i_0 = 0; i_0 < 23; i_0++) {
                    torqueBalancingYoga_B.qj_des[i_0] =
                        torqueBalancingYoga_P.Sm.joints_leftYogaRef[(1 + i_0) * 26 + 25];
                }

                if (rtb_Clock1 > ((torqueBalancingYoga_P.Sm.joints_leftSecondYogaRef[25]
                                   + torqueBalancingYoga_DW.tSwitch)
                                  + torqueBalancingYoga_P.Sm.smoothingTimeSecondYogaLeft)
                                     + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves) {
                    torqueBalancingYoga_DW.state = 5.0;
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                    torqueBalancingYoga_DW.secondYoga = false;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 5.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 0.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            for (i_0 = 0; i_0 < 6; i_0++) {
                LDotDes_0[i_0] = torqueBalancingYoga_B.jointAngles[17 + i_0]
                                 - torqueBalancingYoga_B.qj_des[17 + i_0];
            }

            if (torqueBalancingYoga_norm(LDotDes_0) * 180.0 / 3.1415926535897931
                < torqueBalancingYoga_P.Sm.joints_thresholdNotInContact) {
                for (i_0 = 0; i_0 < 6; i_0++) {
                    LDotDes_0[i_0] = torqueBalancingYoga_B.jointAngles[11 + i_0]
                                     - torqueBalancingYoga_B.qj_des[11 + i_0];
                }

                if (torqueBalancingYoga_norm(LDotDes_0) * 180.0 / 3.1415926535897931
                    < torqueBalancingYoga_P.Sm.joints_thresholdInContact) {
                    torqueBalancingYoga_DW.state = 6.0;
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 6.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 0.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            if (torqueBalancingYoga_B.wR_WBDT[2]
                > torqueBalancingYoga_P.Sm.wrench_thresholdContactOn) {
                torqueBalancingYoga_DW.state = 7.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
            }
        }

        if (torqueBalancingYoga_DW.state == 7.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[0];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[1];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                rtb_impedances[AR_tmp] =
                    torqueBalancingYoga_P.Gain.impedances[(13 * AR_tmp + i_0) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            tmp_1[0] = rtb_CoM6DCoMXYZ2_i[0] - torqueBalancingYoga_B.CoM_des[0];
            tmp_1[1] = rtb_CoM6DCoMXYZ2_i[1] - torqueBalancingYoga_B.CoM_des[1];
            if ((torqueBalancingYoga_norm_i(tmp_1) < 10.0 * torqueBalancingYoga_P.Sm.CoM_threshold)
                && (torqueBalancingYoga_P.Sm.yogaAlsoOnRightFoot
                    && (rtb_Clock1
                        > torqueBalancingYoga_DW.tSwitch + torqueBalancingYoga_P.Sm.tBalancing))) {
                for (i_0 = 0; i_0 < 4; i_0++) {
                    for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] = 0.0;
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6[i_0 << 2]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                    }
                }

                for (i_0 = 0; i_0 < 4; i_0++) {
                    torqueBalancingYoga_DW.w_H_fixedLink[i_0 << 2] = rtb_Switch_jf[i_0 << 2];
                    torqueBalancingYoga_DW.w_H_fixedLink[1 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 1];
                    torqueBalancingYoga_DW.w_H_fixedLink[2 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 2];
                    torqueBalancingYoga_DW.w_H_fixedLink[3 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 3];
                }

                torqueBalancingYoga_mrdivide_o(torqueBalancingYoga_DW.w_H_fixedLink,
                                               torqueBalancingYoga_B.Switch6_g);
                torqueBalancingYoga_DW.state = 8.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
            }
        }

        if (torqueBalancingYoga_DW.state == 8.0) {
            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            fixed_link_CoMDes[0] = torqueBalancingYoga_B.CoM_des[0];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            fixed_link_CoMDes[1] = torqueBalancingYoga_B.CoM_des[1];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            fixed_link_CoMDes[2] = torqueBalancingYoga_B.CoM_des[2];
            fixed_link_CoMDes[3] = 1.0;
            torqueBalancingYoga_mldivide_e(torqueBalancingYoga_DW.w_H_fixedLink, fixed_link_CoMDes);
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            i = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[AR_tmp - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[i - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[AR_tmp + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[i + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[AR_tmp + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[i + 25];
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + i_0) - 1];
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + AR_tmp) - 1];
            }

            if ((std::abs(fixed_link_CoMDes[1] - torqueBalancingYoga_B.SFunction_dn[13])
                 < torqueBalancingYoga_P.Sm.CoM_threshold)
                && (torqueBalancingYoga_B.wL_WBDT[2]
                    < torqueBalancingYoga_P.Sm.wrench_thresholdContactOff)) {
                torqueBalancingYoga_DW.state = 9.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
            }
        }

        if (torqueBalancingYoga_DW.state == 9.0) {
            torqueBalancingYoga_B.constraints[0] = 0.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            if (rtb_Clock1
                > torqueBalancingYoga_DW.tSwitch + torqueBalancingYoga_P.Sm.tBalancingBeforeYoga) {
                torqueBalancingYoga_DW.state = 10.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                if (torqueBalancingYoga_P.Sm.skipYoga) {
                    torqueBalancingYoga_DW.state = 11.0;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 10.0) {
            torqueBalancingYoga_B.constraints[0] = 0.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            for (i = 0; i < 25; i++) {
                if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_rightYogaRef[i]
                                      + torqueBalancingYoga_DW.tSwitch)
                    && (rtb_Clock1 <= torqueBalancingYoga_P.Sm.joints_rightYogaRef[i + 1]
                                          + torqueBalancingYoga_DW.tSwitch)
                    && (!torqueBalancingYoga_DW.secondYoga)) {
                    for (i_0 = 0; i_0 < 23; i_0++) {
                        torqueBalancingYoga_B.qj_des[i_0] =
                            torqueBalancingYoga_P.Sm.joints_rightYogaRef[(1 + i_0) * 26 + i];
                    }
                }
                else {
                    if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_rightSecondYogaRef[i]
                                          + torqueBalancingYoga_DW.tSwitch)
                        && (rtb_Clock1 <= torqueBalancingYoga_P.Sm.joints_rightSecondYogaRef[i + 1]
                                              + torqueBalancingYoga_DW.tSwitch)
                        && torqueBalancingYoga_DW.secondYoga) {
                        for (i_0 = 0; i_0 < 23; i_0++) {
                            torqueBalancingYoga_B.qj_des[i_0] =
                                torqueBalancingYoga_P.Sm.joints_rightYogaRef[(1 + i_0) * 26 + i];
                        }
                    }
                }
            }

            if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_rightYogaRef[25]
                                  + torqueBalancingYoga_DW.tSwitch)
                && (!torqueBalancingYoga_DW.secondYoga)) {
                for (i_0 = 0; i_0 < 23; i_0++) {
                    torqueBalancingYoga_B.qj_des[i_0] =
                        torqueBalancingYoga_P.Sm.joints_rightYogaRef[(1 + i_0) * 26 + 25];
                }

                if ((rtb_Clock1
                     > ((torqueBalancingYoga_P.Sm.joints_rightYogaRef[25]
                         + torqueBalancingYoga_DW.tSwitch)
                        + torqueBalancingYoga_P.Sm
                              .smoothingTimeCoM_Joints[(int32_T) torqueBalancingYoga_DW.state - 1])
                           + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves)
                    && torqueBalancingYoga_P.Sm.repeatYogaMoveset) {
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                    torqueBalancingYoga_DW.secondYoga = true;
                }
                else {
                    if ((rtb_Clock1 > ((torqueBalancingYoga_P.Sm.joints_rightYogaRef[25]
                                        + torqueBalancingYoga_DW.tSwitch)
                                       + torqueBalancingYoga_P.Sm.smoothingTimeCoM_Joints
                                             [(int32_T) torqueBalancingYoga_DW.state - 1])
                                          + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves)
                        && (!torqueBalancingYoga_P.Sm.repeatYogaMoveset)) {
                        torqueBalancingYoga_DW.state = 11.0;
                        torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                        torqueBalancingYoga_DW.secondYoga = false;
                    }
                }
            }

            if ((rtb_Clock1 > torqueBalancingYoga_P.Sm.joints_rightSecondYogaRef[25]
                                  + torqueBalancingYoga_DW.tSwitch)
                && torqueBalancingYoga_DW.secondYoga) {
                for (i_0 = 0; i_0 < 23; i_0++) {
                    torqueBalancingYoga_B.qj_des[i_0] =
                        torqueBalancingYoga_P.Sm.joints_rightYogaRef[(1 + i_0) * 26 + 25];
                }

                if (rtb_Clock1 > ((torqueBalancingYoga_P.Sm.joints_rightSecondYogaRef[25]
                                   + torqueBalancingYoga_DW.tSwitch)
                                  + torqueBalancingYoga_P.Sm.smoothingTimeSecondYogaRight)
                                     + torqueBalancingYoga_P.Sm.joints_pauseBetweenYogaMoves) {
                    torqueBalancingYoga_DW.state = 11.0;
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                    torqueBalancingYoga_DW.secondYoga = false;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 11.0) {
            torqueBalancingYoga_B.constraints[0] = 0.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            for (i_0 = 0; i_0 < 6; i_0++) {
                LDotDes_0[i_0] = torqueBalancingYoga_B.jointAngles[17 + i_0]
                                 - torqueBalancingYoga_B.qj_des[17 + i_0];
            }

            if (torqueBalancingYoga_norm(LDotDes_0) * 180.0 / 3.1415926535897931
                < torqueBalancingYoga_P.Sm.joints_thresholdInContact) {
                for (i_0 = 0; i_0 < 6; i_0++) {
                    LDotDes_0[i_0] = torqueBalancingYoga_B.jointAngles[11 + i_0]
                                     - torqueBalancingYoga_B.qj_des[11 + i_0];
                }

                if (torqueBalancingYoga_norm(LDotDes_0) * 180.0 / 3.1415926535897931
                    < torqueBalancingYoga_P.Sm.joints_thresholdNotInContact) {
                    torqueBalancingYoga_DW.state = 12.0;
                    torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
                }
            }
        }

        if (torqueBalancingYoga_DW.state == 12.0) {
            torqueBalancingYoga_B.constraints[0] = 0.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            torqueBalancingYoga_B.CoM_des[0] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 - 1]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[12];
            torqueBalancingYoga_B.CoM_des[1] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 12]
                                               + torqueBalancingYoga_DW.w_H_fixedLink[13];
            torqueBalancingYoga_B.CoM_des[2] = torqueBalancingYoga_P.Sm.CoM_delta[i_0 + 25]
                                               + torqueBalancingYoga_B.sf_MATLABFunction_f.s0[2];
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.qj_des[i] =
                    torqueBalancingYoga_P.Sm.joints_references[(13 * i + i_0) - 1];
                rtb_impedances[i] = torqueBalancingYoga_P.Gain.impedances[(13 * i + AR_tmp) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            if (torqueBalancingYoga_B.wL_WBDT[2]
                > torqueBalancingYoga_P.Sm.wrench_thresholdContactOn) {
                torqueBalancingYoga_DW.state = 13.0;
                torqueBalancingYoga_DW.tSwitch = rtb_Clock1;
            }
        }

        if (torqueBalancingYoga_DW.state == 13.0) {
            for (i_0 = 0; i_0 < 4; i_0++) {
                for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] = 0.0;
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                    torqueBalancingYoga_B.w_H_b[AR_tmp + (i_0 << 2)] +=
                        torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                        * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                }
            }

            torqueBalancingYoga_B.constraints[0] = 1.0;
            torqueBalancingYoga_B.constraints[1] = 1.0;
            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                rtb_impedances[AR_tmp] =
                    torqueBalancingYoga_P.Gain.impedances[(13 * AR_tmp + i_0) - 1];
            }

            i_0 = (int32_T) torqueBalancingYoga_DW.state;
            AR_tmp = (int32_T) torqueBalancingYoga_DW.state;
            rtb_KPCoM_idx_0 = torqueBalancingYoga_P.Gain.KP_COM[i_0 - 1];
            rtb_KDCoM_idx_0 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp - 1];
            rtb_KPCoM_idx_1 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 12];
            rtb_KDCoM_idx_1 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 12];
            rtb_KPCoM_idx_2 = torqueBalancingYoga_P.Gain.KP_COM[i_0 + 25];
            rtb_KDCoM_idx_2 = torqueBalancingYoga_P.Gain.KD_COM[AR_tmp + 25];
            if ((rtb_Clock1 - torqueBalancingYoga_DW.tSwitch > torqueBalancingYoga_P.Sm.tBalancing)
                && torqueBalancingYoga_P.Sm.yogaInLoop) {
                torqueBalancingYoga_DW.state = 2.0;
                for (i_0 = 0; i_0 < 4; i_0++) {
                    for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] = 0.0;
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6_g[i_0 << 2]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 1]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 2]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                        rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                            torqueBalancingYoga_B.Switch6_g[(i_0 << 2) + 3]
                            * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                    }
                }

                for (i_0 = 0; i_0 < 4; i_0++) {
                    torqueBalancingYoga_DW.w_H_fixedLink[i_0 << 2] = rtb_Switch_jf[i_0 << 2];
                    torqueBalancingYoga_DW.w_H_fixedLink[1 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 1];
                    torqueBalancingYoga_DW.w_H_fixedLink[2 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 2];
                    torqueBalancingYoga_DW.w_H_fixedLink[3 + (i_0 << 2)] =
                        rtb_Switch_jf[(i_0 << 2) + 3];
                }

                torqueBalancingYoga_mrdivide_o(torqueBalancingYoga_DW.w_H_fixedLink,
                                               torqueBalancingYoga_B.Switch6);
                if (torqueBalancingYoga_P.Sm.demoStartsOnRightSupport) {
                    torqueBalancingYoga_DW.state = 8.0;
                    for (i_0 = 0; i_0 < 4; i_0++) {
                        for (AR_tmp = 0; AR_tmp < 4; AR_tmp++) {
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] = 0.0;
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[i_0 << 2]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 1]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 4];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 2]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 8];
                            rtb_Switch_jf[AR_tmp + (i_0 << 2)] +=
                                torqueBalancingYoga_B.Switch6[(i_0 << 2) + 3]
                                * torqueBalancingYoga_DW.w_H_fixedLink[AR_tmp + 12];
                        }
                    }

                    for (i_0 = 0; i_0 < 4; i_0++) {
                        torqueBalancingYoga_DW.w_H_fixedLink[i_0 << 2] = rtb_Switch_jf[i_0 << 2];
                        torqueBalancingYoga_DW.w_H_fixedLink[1 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 1];
                        torqueBalancingYoga_DW.w_H_fixedLink[2 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 2];
                        torqueBalancingYoga_DW.w_H_fixedLink[3 + (i_0 << 2)] =
                            rtb_Switch_jf[(i_0 << 2) + 3];
                    }

                    torqueBalancingYoga_mrdivide_o(torqueBalancingYoga_DW.w_H_fixedLink,
                                                   torqueBalancingYoga_B.Switch6_g);
                }
            }
        }

        if (torqueBalancingYoga_DW.secondYoga && (torqueBalancingYoga_DW.state == 4.0)
            && (rtb_Clock1 >= torqueBalancingYoga_P.Sm.joints_leftSecondYogaRef[1]
                                  + torqueBalancingYoga_DW.tSwitch)) {
            torqueBalancingYoga_B.jointsSmoothingTime =
                torqueBalancingYoga_P.Sm.smoothingTimeSecondYogaLeft;
        }
        else if (torqueBalancingYoga_DW.secondYoga && (torqueBalancingYoga_DW.state == 10.0)
                 && (rtb_Clock1 >= torqueBalancingYoga_P.Sm.joints_rightSecondYogaRef[1]
                                       + torqueBalancingYoga_DW.tSwitch)) {
            torqueBalancingYoga_B.jointsSmoothingTime =
                torqueBalancingYoga_P.Sm.smoothingTimeSecondYogaRight;
        }
        else if ((torqueBalancingYoga_DW.state == 4.0) || (torqueBalancingYoga_DW.state == 10.0)) {
            torqueBalancingYoga_B.jointsSmoothingTime =
                torqueBalancingYoga_P.Sm
                    .smoothingTimeCoM_Joints[(int32_T) torqueBalancingYoga_DW.state - 1]
                * torqueBalancingYoga_P.Sm.scaleFactorSmoothingTime;
        }
        else {
            torqueBalancingYoga_B.jointsSmoothingTime =
                torqueBalancingYoga_P.Sm
                    .smoothingTimeCoM_Joints[(int32_T) torqueBalancingYoga_DW.state - 1];
        }

        /* '<S70>:1:4' */
        torqueBalancingYoga_B.currentState = torqueBalancingYoga_DW.state;

        /* End of MATLAB Function: '<S34>/stateMachineYogaFCN' */

        /* S-Function (WBToolbox): '<S76>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S76>/S-Function

        /* S-Function (WBToolbox): '<S77>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S77>/S-Function

        /* S-Function (WBToolbox): '<S75>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S75>/S-Function

        /* MATLAB Function: '<S65>/Compute Base Velocity' incorporates:
         *  Constant: '<S34>/Constant1'
         */
        /* MATLAB Function 'Robot State and References/State Machine Yoga/Compute State
         * Velocity/Compute Base Velocity': '<S73>:1' */
        /* '<S73>:1:3' */
        for (i_0 = 0; i_0 < 29; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                Gamma_tmp = 6 * i_0 + AR_tmp;
                i = AR_tmp + 12 * i_0;
                Jc[i] = torqueBalancingYoga_B.SFunction_e0[Gamma_tmp]
                        * torqueBalancingYoga_B.constraints[0];
                Jc[i + 6] = torqueBalancingYoga_B.SFunction_da[Gamma_tmp]
                            * torqueBalancingYoga_B.constraints[1];
            }
        }

        for (i_0 = 0; i_0 < 12; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_0[AR_tmp + 6 * i_0] = Jc[12 * AR_tmp + i_0];
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                acc_CoM_des_idx_1 = 0.0;
                for (i = 0; i < 12; i++) {
                    acc_CoM_des_idx_1 += Jc[12 * i_0 + i] * Jc[12 * AR_tmp + i];
                }

                pinvJb[i_0 + 6 * AR_tmp] =
                    (real_T) b_1[6 * AR_tmp + i_0] * torqueBalancingYoga_P.Reg.pinvDamp_nu_b
                    + acc_CoM_des_idx_1;
            }
        }

        torqueBalancingYoga_mldivide(pinvJb, pinvJb_0);
        for (i_0 = 0; i_0 < 12; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_1[AR_tmp + 6 * i_0] = -pinvJb_0[6 * i_0 + AR_tmp];
            }
        }

        for (i = 0; i < 6; i++) {
            torqueBalancingYoga_B.nu_b[i] = 0.0;
            for (i_0 = 0; i_0 < 23; i_0++) {
                Gamma_tmp = i + 6 * i_0;
                pinvJb_2[Gamma_tmp] = 0.0;
                for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                    pinvJb_2[Gamma_tmp] = Jc[(6 + i_0) * 12 + AR_tmp] * pinvJb_1[6 * AR_tmp + i]
                                          + pinvJb_2[6 * i_0 + i];
                }

                torqueBalancingYoga_B.nu_b[i] +=
                    pinvJb_2[6 * i_0 + i] * torqueBalancingYoga_B.SFunction_nr[i_0];
            }

            torqueBalancingYoga_B.Constant1[i] = torqueBalancingYoga_P.Constant1_Value[i];
        }

        /* End of MATLAB Function: '<S65>/Compute Base Velocity' */

        /* SignalConversion: '<S34>/TmpSignal ConversionAtMinimum Jerk Trajectory GeneratorInport1'
         */
        memcpy(&torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[0],
               &rtb_impedances[0],
               23U * sizeof(real_T));
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[23] = rtb_KPCoM_idx_0;
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[26] = rtb_KDCoM_idx_0;
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[24] = rtb_KPCoM_idx_1;
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[27] = rtb_KDCoM_idx_1;
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[25] = rtb_KPCoM_idx_2;
        torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[28] = rtb_KDCoM_idx_2;

        /* S-Function (WBToolbox): '<S34>/Minimum Jerk Trajectory Generator' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S34>/Minimum Jerk Trajectory Generator

        /* Reshape: '<S34>/Reshape1' */
        memcpy(&torqueBalancingYoga_B.Reshape1[0],
               &torqueBalancingYoga_B.w_H_b[0],
               sizeof(real_T) << 4U);
    }

    /* End of Outputs for SubSystem: '<S4>/State Machine Yoga' */

    /* Constant: '<S4>/Constant2' */
    torqueBalancingYoga_B.Constant2 = torqueBalancingYoga_P.Constant2_Value_o;

    /* Outputs for Enabled SubSystem: '<S4>/Internal Coordinator' incorporates:
     *  EnablePort: '<S32>/Enable'
     */
    if (torqueBalancingYoga_B.Compare) {
        if (!torqueBalancingYoga_DW.InternalCoordinator_MODE) {
            torqueBalancingYoga_DW.InternalCoordinator_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.InternalCoordinator_MODE) {
            torqueBalancingYoga_DW.InternalCoordinator_MODE = false;
        }
    }

    if (torqueBalancingYoga_DW.InternalCoordinator_MODE) {
        /* Inport: '<S32>/jointAngles' */
        memcpy(&torqueBalancingYoga_B.jointAngles_p[0],
               &torqueBalancingYoga_B.SFunction_d[0],
               23U * sizeof(real_T));

        /* S-Function (WBToolbox): '<S44>/S-Function' incorporates:
         *  Constant: '<S35>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S44>/S-Function

        /* S-Function (WBToolbox): '<S45>/S-Function' incorporates:
         *  Constant: '<S35>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S45>/S-Function

        /* S-Function (WBToolbox): '<S42>/S-Function' incorporates:
         *  Constant: '<S35>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S42>/S-Function

        /* S-Function (WBToolbox): '<S43>/S-Function' incorporates:
         *  Constant: '<S35>/Constant7'
         */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S43>/S-Function

        /* Switch: '<S35>/Switch' incorporates:
         *  Constant: '<S35>/Constant1'
         */
        if (torqueBalancingYoga_P.Config.LEFT_FOOT_IN_CONTACT_AT_0) {
            memcpy(&rtb_Switch_jf[0], &torqueBalancingYoga_B.SFunction_k0[0], sizeof(real_T) << 4U);
        }
        else {
            memcpy(
                &rtb_Switch_jf[0], &torqueBalancingYoga_B.SFunction_bdr[0], sizeof(real_T) << 4U);
        }

        /* End of Switch: '<S35>/Switch' */

        /* Product: '<S41>/Inv\+' */
        rt_mldivided4x4(torqueBalancingYoga_B.SFunction_jo, rtb_Switch_jf, rtb_imu_H_link_a);

        /* Product: '<S41>/Inv\*  ' */
        rt_mldivided4x4(rtb_Switch_jf, torqueBalancingYoga_B.SFunction_l, rtb_link_H_root_l);

        /* S-Function (WBToolbox): '<S41>/Neck Position' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S41>/Neck Position

        /* RelationalOperator: '<S50>/Compare' incorporates:
         *  Clock: '<S47>/Clock'
         *  Constant: '<S50>/Constant'
         */
        rtb_Compare_g = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const);

        /* MATLAB Function: '<S47>/MATLAB Function' */
        torqueBalancingY_MATLABFunction(rtb_imu_H_link_a,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_os,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_os);

        /* S-Function (WBToolbox): '<S32>/IMU measurements' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S32>/IMU measurements

        /* RelationalOperator: '<S52>/Compare' incorporates:
         *  Clock: '<S48>/Clock'
         *  Constant: '<S52>/Constant'
         */
        rtb_Compare_p = ((&torqueBalancingYoga_M)->Timing.t[0]
                         == torqueBalancingYoga_P.CompareToConstant_const_j);

        /* MATLAB Function: '<S48>/MATLAB Function' */
        torqueBalancin_MATLABFunction_l(torqueBalancingYoga_B.IMUmeasurements,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_l5,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_l5);

        /* Switch: '<S49>/Switch' incorporates:
         *  Constant: '<S49>/Constant'
         *  Constant: '<S49>/USE_IMU4EST_BASE1'
         *  Gain: '<S49>/Gain'
         */
        if (torqueBalancingYoga_P.Config.CORRECT_NECK_IMU) {
            rtb_Switch_j[0] =
                torqueBalancingYoga_P.Gain_Gain * torqueBalancingYoga_B.NeckPosition_c[0];
            rtb_Switch_j[1] =
                torqueBalancingYoga_P.Gain_Gain * torqueBalancingYoga_B.NeckPosition_c[1];
            rtb_Switch_j[2] =
                torqueBalancingYoga_P.Gain_Gain * torqueBalancingYoga_B.NeckPosition_c[2];
        }
        else {
            rtb_Switch_j[0] = torqueBalancingYoga_P.Constant_Value[0];
            rtb_Switch_j[1] = torqueBalancingYoga_P.Constant_Value[1];
            rtb_Switch_j[2] = torqueBalancingYoga_P.Constant_Value[2];
        }

        /* End of Switch: '<S49>/Switch' */

        /* MATLAB Function: '<S41>/fromImuToHomogeousTransformFCN' */
        fromImuToHomogeousTransformFCN(rtb_imu_H_link_a,
                                       torqueBalancingYoga_B.sf_MATLABFunction_os.s0,
                                       rtb_link_H_root_l,
                                       torqueBalancingYoga_B.sf_MATLABFunction_l5.s0,
                                       torqueBalancingYoga_B.IMUmeasurements,
                                       rtb_Switch_j,
                                       &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_f,
                                       &torqueBalancingYoga_P);

        /* Switch: '<S41>/Switch6' incorporates:
         *  Constant: '<S41>/USE_IMU4EST_BASE'
         */
        if (torqueBalancingYoga_P.Config.USE_IMU4EST_BASE) {
            memcpy(&torqueBalancingYoga_B.Switch6_b[0],
                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_f.w_H_b[0],
                   sizeof(real_T) << 4U);
        }
        else {
            memcpy(
                &torqueBalancingYoga_B.Switch6_b[0], &rtb_link_H_root_l[0], sizeof(real_T) << 4U);
        }

        /* End of Switch: '<S41>/Switch6' */

        /* Clock: '<S32>/Clock1' */
        rtb_Clock1 = (&torqueBalancingYoga_M)->Timing.t[0];

        /* S-Function (WBToolbox): '<S57>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S57>/S-Function

        /* S-Function (WBToolbox): '<S58>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S58>/S-Function

        /* Constant: '<S32>/Constant1' */
        torqueBalancingYoga_B.Constant1_c[0] =
            torqueBalancingYoga_P.Config.LEFT_RIGHT_FOOT_IN_CONTACT[0];
        torqueBalancingYoga_B.Constant1_c[1] =
            torqueBalancingYoga_P.Config.LEFT_RIGHT_FOOT_IN_CONTACT[1];

        /* S-Function (WBToolbox): '<S56>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S56>/S-Function

        /* MATLAB Function: '<S36>/Compute Base Velocity' */
        /* MATLAB Function 'Robot State and References/Internal Coordinator/Compute State
         * Velocity/Compute Base Velocity': '<S54>:1' */
        /* '<S54>:1:3' */
        for (i_0 = 0; i_0 < 29; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                Gamma_tmp = 6 * i_0 + AR_tmp;
                i = AR_tmp + 12 * i_0;
                Jc[i] = torqueBalancingYoga_B.SFunction_p[Gamma_tmp]
                        * torqueBalancingYoga_B.Constant1_c[0];
                Jc[i + 6] = torqueBalancingYoga_B.SFunction_i0[Gamma_tmp]
                            * torqueBalancingYoga_B.Constant1_c[1];
            }
        }

        for (i_0 = 0; i_0 < 12; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_0[AR_tmp + 6 * i_0] = Jc[12 * AR_tmp + i_0];
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                acc_CoM_des_idx_1 = 0.0;
                for (i = 0; i < 12; i++) {
                    acc_CoM_des_idx_1 += Jc[12 * i_0 + i] * Jc[12 * AR_tmp + i];
                }

                pinvJb[i_0 + 6 * AR_tmp] =
                    (real_T) b_2[6 * AR_tmp + i_0] * torqueBalancingYoga_P.Reg.pinvDamp_nu_b
                    + acc_CoM_des_idx_1;
            }
        }

        torqueBalancingYoga_mldivide(pinvJb, pinvJb_0);
        for (i_0 = 0; i_0 < 12; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_1[AR_tmp + 6 * i_0] = -pinvJb_0[6 * i_0 + AR_tmp];
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            torqueBalancingYoga_B.nu_b_e[i_0] = 0.0;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                Gamma_tmp = i_0 + 6 * AR_tmp;
                pinvJb_2[Gamma_tmp] = 0.0;
                for (i = 0; i < 12; i++) {
                    pinvJb_2[Gamma_tmp] = Jc[(6 + AR_tmp) * 12 + i] * pinvJb_1[6 * i + i_0]
                                          + pinvJb_2[6 * AR_tmp + i_0];
                }

                torqueBalancingYoga_B.nu_b_e[i_0] +=
                    pinvJb_2[6 * AR_tmp + i_0] * torqueBalancingYoga_B.SFunction_m[AR_tmp];
            }
        }

        /* End of MATLAB Function: '<S36>/Compute Base Velocity' */

        /* S-Function (WBToolbox): '<S63>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S63>/S-Function

        /* Selector: '<S40>/CoM6D ->  CoMXYZ2' */
        rtb_CoM6DCoMXYZ2_b[0] = torqueBalancingYoga_B.SFunction_ak[12];
        rtb_CoM6DCoMXYZ2_b[1] = torqueBalancingYoga_B.SFunction_ak[13];
        rtb_CoM6DCoMXYZ2_b[2] = torqueBalancingYoga_B.SFunction_ak[14];

        /* RelationalOperator: '<S61>/Compare' incorporates:
         *  Clock: '<S39>/Clock'
         *  Constant: '<S61>/Constant'
         */
        rtb_Compare_f4 = ((&torqueBalancingYoga_M)->Timing.t[0]
                          == torqueBalancingYoga_P.CompareToConstant_const_o);

        /* MATLAB Function: '<S39>/MATLAB Function' */
        torqueBalancin_MATLABFunction_a(rtb_CoM6DCoMXYZ2_b,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_a,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_a);

        /* MATLAB Function: '<S32>/Reference Generator CoM' */
        /* MATLAB Function 'Robot State and References/Internal Coordinator/Reference Generator
         * CoM': '<S37>:1' */
        /* '<S37>:1:3' */
        rtb_KPCoM_idx_0 = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[0];
        rtb_KDCoM_idx_0 = 0.0;
        acc_CoM_des_idx_0 = 0.0;
        rtb_KPCoM_idx_1 = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[1];
        rtb_KDCoM_idx_1 = 0.0;
        acc_CoM_des_idx_1 = 0.0;
        rtb_KPCoM_idx_2 = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[2];
        rtb_KDCoM_idx_2 = 0.0;
        acc_CoM_des_idx_2 = 0.0;
        if (torqueBalancingYoga_P.Config.amplitudeOfOscillation != 0.0) {
            if (rtb_Clock1 > torqueBalancingYoga_P.Config.noOscillationTime) {
                acc_CoM_des_idx_0 = torqueBalancingYoga_P.Config.amplitudeOfOscillation;
            }
            else {
                acc_CoM_des_idx_0 = 0.0;
            }

            acc_CoM_des_idx_1 =
                6.2831853071795862 * torqueBalancingYoga_P.Config.frequencyOfOscillation;
            rtb_KDCoM_idx_2 = std::sin(acc_CoM_des_idx_1 * rtb_Clock1) * acc_CoM_des_idx_0;
            rtb_KPCoM_idx_0 =
                rtb_KDCoM_idx_2 * torqueBalancingYoga_P.Config.directionOfOscillation[0]
                + torqueBalancingYoga_B.sf_MATLABFunction_a.s0[0];
            rtb_KPCoM_idx_1 =
                rtb_KDCoM_idx_2 * torqueBalancingYoga_P.Config.directionOfOscillation[1]
                + torqueBalancingYoga_B.sf_MATLABFunction_a.s0[1];
            rtb_KPCoM_idx_2 =
                rtb_KDCoM_idx_2 * torqueBalancingYoga_P.Config.directionOfOscillation[2]
                + torqueBalancingYoga_B.sf_MATLABFunction_a.s0[2];
            rtb_KDCoM_idx_2 = acc_CoM_des_idx_0 * 2.0 * 3.1415926535897931
                              * torqueBalancingYoga_P.Config.frequencyOfOscillation
                              * std::cos(acc_CoM_des_idx_1 * rtb_Clock1);
            rtb_KDCoM_idx_0 =
                rtb_KDCoM_idx_2 * torqueBalancingYoga_P.Config.directionOfOscillation[0];
            rtb_KDCoM_idx_1 =
                rtb_KDCoM_idx_2 * torqueBalancingYoga_P.Config.directionOfOscillation[1];
            rtb_KDCoM_idx_2 *= torqueBalancingYoga_P.Config.directionOfOscillation[2];
            rtb_Clock1 = acc_CoM_des_idx_1 * acc_CoM_des_idx_1 * -acc_CoM_des_idx_0
                         * std::sin(acc_CoM_des_idx_1 * rtb_Clock1);
            acc_CoM_des_idx_0 = rtb_Clock1 * torqueBalancingYoga_P.Config.directionOfOscillation[0];
            acc_CoM_des_idx_1 = rtb_Clock1 * torqueBalancingYoga_P.Config.directionOfOscillation[1];
            acc_CoM_des_idx_2 = rtb_Clock1 * torqueBalancingYoga_P.Config.directionOfOscillation[2];
        }

        /* Reshape: '<S32>/Reshape1' */
        memcpy(&torqueBalancingYoga_B.Reshape1_p[0],
               &torqueBalancingYoga_B.Switch6_b[0],
               sizeof(real_T) << 4U);

        /* Switch: '<S32>/Switch1' incorporates:
         *  Constant: '<S32>/Constant'
         *  Constant: '<S32>/Constant2'
         *  MATLAB Function: '<S32>/Reference Generator CoM'
         */
        if (torqueBalancingYoga_P.Config.DEMO_MOVEMENTS) {
            torqueBalancingYoga_B.Switch1[0] = rtb_KPCoM_idx_0;
            torqueBalancingYoga_B.Switch1[3] = rtb_KDCoM_idx_0;
            torqueBalancingYoga_B.Switch1[6] = acc_CoM_des_idx_0;
            torqueBalancingYoga_B.Switch1[1] = rtb_KPCoM_idx_1;
            torqueBalancingYoga_B.Switch1[4] = rtb_KDCoM_idx_1;
            torqueBalancingYoga_B.Switch1[7] = acc_CoM_des_idx_1;
            torqueBalancingYoga_B.Switch1[2] = rtb_KPCoM_idx_2;
            torqueBalancingYoga_B.Switch1[5] = rtb_KDCoM_idx_2;
            torqueBalancingYoga_B.Switch1[8] = acc_CoM_des_idx_2;
        }
        else {
            torqueBalancingYoga_B.Switch1[0] = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[0];
            torqueBalancingYoga_B.Switch1[3] = torqueBalancingYoga_P.Constant2_Value[0];
            torqueBalancingYoga_B.Switch1[6] = torqueBalancingYoga_P.Constant2_Value[3];
            torqueBalancingYoga_B.Switch1[1] = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[1];
            torqueBalancingYoga_B.Switch1[4] = torqueBalancingYoga_P.Constant2_Value[1];
            torqueBalancingYoga_B.Switch1[7] = torqueBalancingYoga_P.Constant2_Value[4];
            torqueBalancingYoga_B.Switch1[2] = torqueBalancingYoga_B.sf_MATLABFunction_a.s0[2];
            torqueBalancingYoga_B.Switch1[5] = torqueBalancingYoga_P.Constant2_Value[2];
            torqueBalancingYoga_B.Switch1[8] = torqueBalancingYoga_P.Constant2_Value[5];
        }

        /* End of Switch: '<S32>/Switch1' */

        /* RelationalOperator: '<S59>/Compare' incorporates:
         *  Clock: '<S38>/Clock'
         *  Constant: '<S59>/Constant'
         */
        rtb_Compare_fn = ((&torqueBalancingYoga_M)->Timing.t[0]
                          == torqueBalancingYoga_P.CompareToConstant_const_ok);

        /* MATLAB Function: '<S38>/MATLAB Function' */
        torqueBalancin_MATLABFunction_k(torqueBalancingYoga_B.jointAngles_p,
                                        &torqueBalancingYoga_B.sf_MATLABFunction_k,
                                        &torqueBalancingYoga_DW.sf_MATLABFunction_k);

        /* Constant: '<S32>/Constant3' */
        memcpy(&torqueBalancingYoga_B.Constant3[0],
               &torqueBalancingYoga_P.Constant3_Value[0],
               23U * sizeof(real_T));

        /* Constant: '<S32>/Constant4' */
        torqueBalancingYoga_B.Constant4 = torqueBalancingYoga_P.Constant4_Value;

        /* Constant: '<S32>/joints.smoothingTime' */
        torqueBalancingYoga_B.jointssmoothingTime =
            torqueBalancingYoga_P.Sm.smoothingTimeCoM_Joints[0];

        /* Constant: '<S32>/Constant5' */
        torqueBalancingYoga_B.Constant5[0] = torqueBalancingYoga_P.Constant5_Value[0];

        /* Constant: '<S32>/Constant6' */
        torqueBalancingYoga_B.Constant6[0] = torqueBalancingYoga_P.Constant6_Value[0];

        /* Constant: '<S32>/Constant5' */
        torqueBalancingYoga_B.Constant5[1] = torqueBalancingYoga_P.Constant5_Value[1];

        /* Constant: '<S32>/Constant6' */
        torqueBalancingYoga_B.Constant6[1] = torqueBalancingYoga_P.Constant6_Value[1];

        /* Constant: '<S32>/Constant5' */
        torqueBalancingYoga_B.Constant5[2] = torqueBalancingYoga_P.Constant5_Value[2];

        /* Constant: '<S32>/Constant6' */
        torqueBalancingYoga_B.Constant6[2] = torqueBalancingYoga_P.Constant6_Value[2];
    }

    /* End of Outputs for SubSystem: '<S4>/Internal Coordinator' */

    /* MultiPortSwitch: '<S4>/Multiport Switch1' */
    if (torqueBalancingYoga_B.Constant2 == 1) {
        memcpy(&torqueBalancingYoga_B.MultiportSwitch1[0],
               &torqueBalancingYoga_B.Switch1[0],
               9U * sizeof(real_T));
        torqueBalancingYoga_B.MultiportSwitch1[32] = torqueBalancingYoga_B.Constant1_c[0];
        torqueBalancingYoga_B.MultiportSwitch1[33] = torqueBalancingYoga_B.Constant1_c[1];
        for (i = 0; i < 23; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 9] =
                torqueBalancingYoga_B.sf_MATLABFunction_k.s0[i];
            torqueBalancingYoga_B.MultiportSwitch1[i + 34] = torqueBalancingYoga_B.Constant3[i];
        }

        torqueBalancingYoga_B.MultiportSwitch1[57] = torqueBalancingYoga_B.Constant4;
        for (i = 0; i < 6; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 58] = torqueBalancingYoga_B.nu_b_e[i];
        }

        for (i = 0; i < 23; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 64] = torqueBalancingYoga_B.SFunction_m[i];
            torqueBalancingYoga_B.MultiportSwitch1[i + 87] = torqueBalancingYoga_B.jointAngles_p[i];
        }

        torqueBalancingYoga_B.MultiportSwitch1[110] = torqueBalancingYoga_B.jointssmoothingTime;
        torqueBalancingYoga_B.MultiportSwitch1[111] = torqueBalancingYoga_B.Constant5[0];
        torqueBalancingYoga_B.MultiportSwitch1[114] = torqueBalancingYoga_B.Constant6[0];
        torqueBalancingYoga_B.MultiportSwitch1[112] = torqueBalancingYoga_B.Constant5[1];
        torqueBalancingYoga_B.MultiportSwitch1[115] = torqueBalancingYoga_B.Constant6[1];
        torqueBalancingYoga_B.MultiportSwitch1[113] = torqueBalancingYoga_B.Constant5[2];
        torqueBalancingYoga_B.MultiportSwitch1[116] = torqueBalancingYoga_B.Constant6[2];
        memcpy(&torqueBalancingYoga_B.MultiportSwitch1[117],
               &torqueBalancingYoga_B.Reshape1_p[0],
               sizeof(real_T) << 4U);
    }
    else {
        torqueBalancingYoga_B.MultiportSwitch1[0] = torqueBalancingYoga_B.CoM_des[0];
        torqueBalancingYoga_B.MultiportSwitch1[1] = torqueBalancingYoga_B.CoM_des[1];
        torqueBalancingYoga_B.MultiportSwitch1[2] = torqueBalancingYoga_B.CoM_des[2];
        for (i = 0; i < 6; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 3] = torqueBalancingYoga_B.Constant1[i];
        }

        torqueBalancingYoga_B.MultiportSwitch1[32] = torqueBalancingYoga_B.constraints[0];
        torqueBalancingYoga_B.MultiportSwitch1[33] = torqueBalancingYoga_B.constraints[1];
        for (i = 0; i < 23; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 9] = torqueBalancingYoga_B.qj_des[i];
            torqueBalancingYoga_B.MultiportSwitch1[i + 34] =
                torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator[i];
        }

        torqueBalancingYoga_B.MultiportSwitch1[57] = torqueBalancingYoga_B.currentState;
        for (i = 0; i < 6; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 58] = torqueBalancingYoga_B.nu_b[i];
        }

        for (i = 0; i < 23; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 64] = torqueBalancingYoga_B.SFunction_nr[i];
            torqueBalancingYoga_B.MultiportSwitch1[i + 87] = torqueBalancingYoga_B.jointAngles[i];
        }

        torqueBalancingYoga_B.MultiportSwitch1[110] = torqueBalancingYoga_B.jointsSmoothingTime;
        for (i = 0; i < 6; i++) {
            torqueBalancingYoga_B.MultiportSwitch1[i + 111] =
                torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator[i + 23];
        }

        memcpy(&torqueBalancingYoga_B.MultiportSwitch1[117],
               &torqueBalancingYoga_B.Reshape1[0],
               sizeof(real_T) << 4U);
    }

    /* End of MultiPortSwitch: '<S4>/Multiport Switch1' */

    /* S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator1' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator1

    /* Switch: '<S33>/Switch5' incorporates:
     *  Constant: '<S33>/SMOOTH_DES_COM2'
     */
    if (torqueBalancingYoga_P.Config.SMOOTH_JOINT_DES) {
        memcpy(&torqueBalancingYoga_B.Switch5[0],
               &torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator1[0],
               23U * sizeof(real_T));
    }
    else {
        memcpy(&torqueBalancingYoga_B.Switch5[0],
               &torqueBalancingYoga_B.MultiportSwitch1[9],
               23U * sizeof(real_T));
    }

    /* End of Switch: '<S33>/Switch5' */

    /* S-Function (WBToolbox): '<S18>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S18>/S-Function

    /* MATLAB Function: '<S15>/Add motor reflected inertias' */
    /* MATLAB Function 'Dynamics and Kinematics/Dynamics/Add motors reflected inertias/Add motor
     * reflected inertias': '<S19>:1' */
    /* '<S19>:1:3' */
    if (torqueBalancingYoga_P.Config.USE_MOTOR_REFLECTED_INERTIA) {
        torqueBala_computeMotorsInertia(&torqueBalancingYoga_P.Config,
                                        torqueBalancingYoga_B.invTGamma);
        for (i_0 = 0; i_0 < 29; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                torqueBalancingYoga_B.M_with_inertia[AR_tmp + 29 * i_0] = 0.0;
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            memset(&torqueBalancingYoga_B.M_with_inertia[i_0 * 29 + 6], 0, 23U * sizeof(real_T));
        }

        for (i_0 = 0; i_0 < 23; i_0++) {
            memcpy(&torqueBalancingYoga_B.M_with_inertia[i_0 * 29 + 180],
                   &torqueBalancingYoga_B.invTGamma[i_0 * 23],
                   23U * sizeof(real_T));
        }
    }
    else {
        memset(&torqueBalancingYoga_B.M_with_inertia[0], 0, 841U * sizeof(real_T));
    }

    for (i_0 = 0; i_0 < 841; i_0++) {
        torqueBalancingYoga_B.M_with_inertia[i_0] += torqueBalancingYoga_B.SFunction_b[i_0];
    }

    /* End of MATLAB Function: '<S15>/Add motor reflected inertias' */

    /* Gain: '<S17>/Gain' */
    for (i = 0; i < 23; i++) {
        torqueBalancingYoga_B.Gain[i] =
            torqueBalancingYoga_B.MultiportSwitch1[i + 64] * torqueBalancingYoga_P.Gain_Gain_o;
    }

    /* End of Gain: '<S17>/Gain' */

    /* S-Function (WBToolbox): '<S17>/S-Function' incorporates:
     *  Constant: '<S17>/Constant'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S17>/S-Function

    /* S-Function (WBToolbox): '<S16>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S16>/S-Function

    /* S-Function (WBToolbox): '<S122>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S122>/S-Function

    /* S-Function (WBToolbox): '<S119>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S119>/S-Function

    /* Product: '<S118>/Inv\+' */
    rt_mldivided4x4(
        torqueBalancingYoga_B.SFunction_bd, torqueBalancingYoga_B.SFunction_br, rtb_imu_H_link);

    /* RelationalOperator: '<S128>/Compare' incorporates:
     *  Clock: '<S125>/Clock'
     *  Constant: '<S128>/Constant'
     */
    rtb_Compare =
        ((&torqueBalancingYoga_M)->Timing.t[0] == torqueBalancingYoga_P.CompareToConstant_const_h);

    /* MATLAB Function: '<S125>/MATLAB Function' */
    torqueBalancingY_MATLABFunction(rtb_imu_H_link,
                                    &torqueBalancingYoga_B.sf_MATLABFunction,
                                    &torqueBalancingYoga_DW.sf_MATLABFunction);

    /* S-Function (WBToolbox): '<S123>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S123>/S-Function

    /* Product: '<S118>/Inv\*  ' */
    rt_mldivided4x4(
        torqueBalancingYoga_B.SFunction_br, torqueBalancingYoga_B.SFunction_d2, rtb_link_H_root);

    /* S-Function (WBToolbox): '<S109>/inertial' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.inertial_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.inertial_PWORK.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S109>/inertial

    /* RelationalOperator: '<S130>/Compare' incorporates:
     *  Clock: '<S126>/Clock'
     *  Constant: '<S130>/Constant'
     */
    rtb_Compare_c =
        ((&torqueBalancingYoga_M)->Timing.t[0] == torqueBalancingYoga_P.CompareToConstant_const_c);

    /* MATLAB Function: '<S126>/MATLAB Function' */
    torqueBalancin_MATLABFunction_l(torqueBalancingYoga_B.inertial,
                                    &torqueBalancingYoga_B.sf_MATLABFunction_h,
                                    &torqueBalancingYoga_DW.sf_MATLABFunction_h);

    /* S-Function (WBToolbox): '<S118>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S118>/Neck Position

    /* Switch: '<S127>/Switch' incorporates:
     *  Constant: '<S127>/Constant'
     *  Constant: '<S127>/USE_IMU4EST_BASE1'
     *  Gain: '<S127>/Gain'
     */
    if (torqueBalancingYoga_P.Config.CORRECT_NECK_IMU) {
        rtb_Switch[0] = torqueBalancingYoga_P.Gain_Gain_hk * torqueBalancingYoga_B.NeckPosition[0];
        rtb_Switch[1] = torqueBalancingYoga_P.Gain_Gain_hk * torqueBalancingYoga_B.NeckPosition[1];
        rtb_Switch[2] = torqueBalancingYoga_P.Gain_Gain_hk * torqueBalancingYoga_B.NeckPosition[2];
    }
    else {
        rtb_Switch[0] = torqueBalancingYoga_P.Constant_Value_o[0];
        rtb_Switch[1] = torqueBalancingYoga_P.Constant_Value_o[1];
        rtb_Switch[2] = torqueBalancingYoga_P.Constant_Value_o[2];
    }

    /* End of Switch: '<S127>/Switch' */

    /* MATLAB Function: '<S118>/fromImuToHomogeousTransformFCN' */
    fromImuToHomogeousTransformFCN(rtb_imu_H_link,
                                   torqueBalancingYoga_B.sf_MATLABFunction.s0,
                                   rtb_link_H_root,
                                   torqueBalancingYoga_B.sf_MATLABFunction_h.s0,
                                   torqueBalancingYoga_B.inertial,
                                   rtb_Switch,
                                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransformF,
                                   &torqueBalancingYoga_P);

    /* MATLAB Function: '<S109>/choose base to world transform' */
    /* MATLAB Function 'controller_QP/Compute angular momentum integral/choose base to world
     * transform': '<S117>:1' */
    /* '<S117>:1:8' */
    rtb_Clock1 = 1.0;
    if ((torqueBalancingYoga_B.MultiportSwitch1[57] == 9.0)
        || (torqueBalancingYoga_B.MultiportSwitch1[57] == 10.0)
        || (torqueBalancingYoga_B.MultiportSwitch1[57] == 11.0)
        || (torqueBalancingYoga_B.MultiportSwitch1[57] == 12.0)) {
        /* '<S117>:1:10' */
        /* '<S117>:1:12' */
        rtb_Clock1 = 0.0;
    }

    /* End of MATLAB Function: '<S109>/choose base to world transform' */

    /* S-Function (WBToolbox): '<S132>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S132>/S-Function

    /* S-Function (WBToolbox): '<S121>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S121>/S-Function

    /* Product: '<S120>/Inv\+' */
    rt_mldivided4x4(
        torqueBalancingYoga_B.SFunction_o, torqueBalancingYoga_B.SFunction_a, rtb_imu_H_link_f);

    /* RelationalOperator: '<S138>/Compare' incorporates:
     *  Clock: '<S135>/Clock'
     *  Constant: '<S138>/Constant'
     */
    rtb_Compare_i =
        ((&torqueBalancingYoga_M)->Timing.t[0] == torqueBalancingYoga_P.CompareToConstant_const_ea);

    /* MATLAB Function: '<S135>/MATLAB Function' */
    torqueBalancingY_MATLABFunction(rtb_imu_H_link_f,
                                    &torqueBalancingYoga_B.sf_MATLABFunction_o,
                                    &torqueBalancingYoga_DW.sf_MATLABFunction_o);

    /* S-Function (WBToolbox): '<S133>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S133>/S-Function

    /* Product: '<S120>/Inv\*  ' */
    rt_mldivided4x4(
        torqueBalancingYoga_B.SFunction_a, torqueBalancingYoga_B.SFunction_c2, rtb_link_H_root_b);

    /* RelationalOperator: '<S140>/Compare' incorporates:
     *  Clock: '<S136>/Clock'
     *  Constant: '<S140>/Constant'
     */
    rtb_Compare_h =
        ((&torqueBalancingYoga_M)->Timing.t[0] == torqueBalancingYoga_P.CompareToConstant_const_e5);

    /* MATLAB Function: '<S136>/MATLAB Function' */
    torqueBalancin_MATLABFunction_l(torqueBalancingYoga_B.inertial,
                                    &torqueBalancingYoga_B.sf_MATLABFunction_l,
                                    &torqueBalancingYoga_DW.sf_MATLABFunction_l);

    /* S-Function (WBToolbox): '<S120>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S120>/Neck Position

    /* Switch: '<S137>/Switch' incorporates:
     *  Constant: '<S137>/Constant'
     *  Constant: '<S137>/USE_IMU4EST_BASE1'
     *  Gain: '<S137>/Gain'
     */
    if (torqueBalancingYoga_P.Config.CORRECT_NECK_IMU) {
        rtb_Switch_l[0] =
            torqueBalancingYoga_P.Gain_Gain_c * torqueBalancingYoga_B.NeckPosition_m[0];
        rtb_Switch_l[1] =
            torqueBalancingYoga_P.Gain_Gain_c * torqueBalancingYoga_B.NeckPosition_m[1];
        rtb_Switch_l[2] =
            torqueBalancingYoga_P.Gain_Gain_c * torqueBalancingYoga_B.NeckPosition_m[2];
    }
    else {
        rtb_Switch_l[0] = torqueBalancingYoga_P.Constant_Value_i[0];
        rtb_Switch_l[1] = torqueBalancingYoga_P.Constant_Value_i[1];
        rtb_Switch_l[2] = torqueBalancingYoga_P.Constant_Value_i[2];
    }

    /* End of Switch: '<S137>/Switch' */

    /* MATLAB Function: '<S120>/fromImuToHomogeousTransformFCN' */
    fromImuToHomogeousTransformFCN(rtb_imu_H_link_f,
                                   torqueBalancingYoga_B.sf_MATLABFunction_o.s0,
                                   rtb_link_H_root_b,
                                   torqueBalancingYoga_B.sf_MATLABFunction_l.s0,
                                   torqueBalancingYoga_B.inertial,
                                   rtb_Switch_l,
                                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_e,
                                   &torqueBalancingYoga_P);

    /* Switch: '<S109>/Switch' incorporates:
     *  Constant: '<S118>/USE_IMU4EST_BASE'
     *  Constant: '<S120>/USE_IMU4EST_BASE'
     *  Switch: '<S118>/Switch6'
     *  Switch: '<S120>/Switch6'
     */
    if (rtb_Clock1 > torqueBalancingYoga_P.Switch_Threshold) {
        if (torqueBalancingYoga_P.Config.USE_IMU4EST_BASE) {
            memcpy(&torqueBalancingYoga_B.Switch[0],
                   &torqueBalancingYoga_B.sf_fromImuToHomogeousTransformF.w_H_b[0],
                   sizeof(real_T) << 4U);
        }
        else {
            memcpy(&torqueBalancingYoga_B.Switch[0], &rtb_link_H_root[0], sizeof(real_T) << 4U);
        }
    }
    else if (torqueBalancingYoga_P.Config.USE_IMU4EST_BASE) {
        /* Switch: '<S120>/Switch6' incorporates:
         *  Switch: '<S118>/Switch6'
         */
        memcpy(&torqueBalancingYoga_B.Switch[0],
               &torqueBalancingYoga_B.sf_fromImuToHomogeousTransfor_e.w_H_b[0],
               sizeof(real_T) << 4U);
    }
    else {
        memcpy(&torqueBalancingYoga_B.Switch[0], &rtb_link_H_root_b[0], sizeof(real_T) << 4U);
    }

    /* End of Switch: '<S109>/Switch' */

    /* Sum: '<S109>/Sum' */
    for (i = 0; i < 23; i++) {
        torqueBalancingYoga_B.Sum[i] =
            torqueBalancingYoga_B.MultiportSwitch1[i + 87] - torqueBalancingYoga_B.Switch5[i];
    }

    /* End of Sum: '<S109>/Sum' */

    /* S-Function (WBToolbox): '<S114>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S114>/S-Function

    /* S-Function (WBToolbox): '<S115>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S115>/S-Function

    /* MATLAB Function: '<S109>/References for L' */
    /* MATLAB Function 'controller_QP/Compute angular momentum integral/References for L':
     * '<S116>:1' */
    /* '<S116>:1:5' */
    if (torqueBalancingYoga_B.MultiportSwitch1[32] == 1.0) {
        /* '<S116>:1:3' */
        /* '<S116>:1:5' */
        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                Gamma_tmp = 6 * AR_tmp + i_0;
                pinvJb[AR_tmp + 6 * i_0] = torqueBalancingYoga_B.SFunction_e[Gamma_tmp];
                acc_CoM_des_idx_1 = 0.0;
                for (i = 0; i < 6; i++) {
                    acc_CoM_des_idx_1 += torqueBalancingYoga_B.SFunction_e[6 * i_0 + i]
                                         * torqueBalancingYoga_B.SFunction_e[6 * AR_tmp + i];
                }

                pinvJb_3[i_0 + 6 * AR_tmp] =
                    (real_T) b_0[Gamma_tmp] * torqueBalancingYoga_P.Reg.pinvDamp_nu_b
                    + acc_CoM_des_idx_1;
            }
        }

        torqueBalancingYoga_mldivide_h(pinvJb_3, pinvJb);

        /* '<S116>:1:6' */
        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_3[AR_tmp + 6 * i_0] = -pinvJb[6 * i_0 + AR_tmp];
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            torqueBalancingYoga_B.nu_b_equivalent[i_0] = 0.0;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                Gamma_tmp = i_0 + 6 * AR_tmp;
                pinvJb_2[Gamma_tmp] = 0.0;
                for (i = 0; i < 6; i++) {
                    pinvJb_2[Gamma_tmp] = torqueBalancingYoga_B.SFunction_e[(6 + AR_tmp) * 6 + i]
                                              * pinvJb_3[6 * i + i_0]
                                          + pinvJb_2[6 * AR_tmp + i_0];
                }

                torqueBalancingYoga_B.nu_b_equivalent[i_0] +=
                    pinvJb_2[6 * AR_tmp + i_0] * torqueBalancingYoga_B.Sum[AR_tmp];
            }
        }
    }
    else {
        /* '<S116>:1:9' */
        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                Gamma_tmp = 6 * AR_tmp + i_0;
                pinvJb[AR_tmp + 6 * i_0] = torqueBalancingYoga_B.SFunction_cf[Gamma_tmp];
                acc_CoM_des_idx_1 = 0.0;
                for (i = 0; i < 6; i++) {
                    acc_CoM_des_idx_1 += torqueBalancingYoga_B.SFunction_cf[6 * i_0 + i]
                                         * torqueBalancingYoga_B.SFunction_cf[6 * AR_tmp + i];
                }

                pinvJb_3[i_0 + 6 * AR_tmp] =
                    (real_T) b_0[Gamma_tmp] * torqueBalancingYoga_P.Reg.pinvDamp_nu_b
                    + acc_CoM_des_idx_1;
            }
        }

        torqueBalancingYoga_mldivide_h(pinvJb_3, pinvJb);

        /* '<S116>:1:10' */
        for (i_0 = 0; i_0 < 6; i_0++) {
            for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                pinvJb_3[AR_tmp + 6 * i_0] = -pinvJb[6 * i_0 + AR_tmp];
            }
        }

        for (i_0 = 0; i_0 < 6; i_0++) {
            torqueBalancingYoga_B.nu_b_equivalent[i_0] = 0.0;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                Gamma_tmp = i_0 + 6 * AR_tmp;
                pinvJb_2[Gamma_tmp] = 0.0;
                for (i = 0; i < 6; i++) {
                    pinvJb_2[Gamma_tmp] = torqueBalancingYoga_B.SFunction_cf[(6 + AR_tmp) * 6 + i]
                                              * pinvJb_3[6 * i + i_0]
                                          + pinvJb_2[6 * AR_tmp + i_0];
                }

                torqueBalancingYoga_B.nu_b_equivalent[i_0] +=
                    pinvJb_2[6 * AR_tmp + i_0] * torqueBalancingYoga_B.Sum[AR_tmp];
            }
        }
    }

    /* End of MATLAB Function: '<S109>/References for L' */

    /* S-Function (WBToolbox): '<S112>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S112>/S-Function

    /* S-Function (WBToolbox): '<S25>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S25>/S-Function

    /* S-Function (WBToolbox): '<S26>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S26>/S-Function

    /* S-Function (WBToolbox): '<S23>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S23>/S-Function

    /* S-Function (WBToolbox): '<S24>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S24>/S-Function

    /* S-Function (WBToolbox): '<S20>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S20>/S-Function

    /* S-Function (WBToolbox): '<S21>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S21>/S-Function

    /* S-Function (WBToolbox): '<S28>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S28>/S-Function

    /* S-Function (WBToolbox): '<S22>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S22>/S-Function

    /* S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator2' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator2

    /* Switch: '<S33>/Switch3' incorporates:
     *  Constant: '<S33>/SMOOTH_DES_COM'
     */
    if (torqueBalancingYoga_P.Config.SMOOTH_COM_DES) {
        rtb_Switch3[0] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator2[0];
        rtb_Switch3[3] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_d[0];
        rtb_Switch3[6] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_b[0];
        rtb_Switch3[1] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator2[1];
        rtb_Switch3[4] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_d[1];
        rtb_Switch3[7] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_b[1];
        rtb_Switch3[2] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator2[2];
        rtb_Switch3[5] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_d[2];
        rtb_Switch3[8] = torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_b[2];
    }
    else {
        rtb_Switch3[0] = torqueBalancingYoga_B.MultiportSwitch1[0];
        rtb_Switch3[3] = torqueBalancingYoga_B.MultiportSwitch1[3];
        rtb_Switch3[6] = torqueBalancingYoga_B.MultiportSwitch1[6];
        rtb_Switch3[1] = torqueBalancingYoga_B.MultiportSwitch1[1];
        rtb_Switch3[4] = torqueBalancingYoga_B.MultiportSwitch1[4];
        rtb_Switch3[7] = torqueBalancingYoga_B.MultiportSwitch1[7];
        rtb_Switch3[2] = torqueBalancingYoga_B.MultiportSwitch1[2];
        rtb_Switch3[5] = torqueBalancingYoga_B.MultiportSwitch1[5];
        rtb_Switch3[8] = torqueBalancingYoga_B.MultiportSwitch1[8];
    }

    /* End of Switch: '<S33>/Switch3' */

    /* MATLAB Function: '<S6>/Balancing Controller ' incorporates:
     *  Constant: '<S6>/Constant'
     *  Constant: '<S6>/Constant1'
     *  Selector: '<S27>/CoM6D ->  CoMXYZ2'
     */
    /* MATLAB Function 'controller_QP/Balancing Controller ': '<S108>:1' */
    /* '<S108>:1:9' */
    gravityWrench[0] = 0.0;
    gravityWrench[1] = 0.0;
    gravityWrench[2] = -torqueBalancingYoga_B.M_with_inertia[0] * 9.81;
    gravityWrench[3] = 0.0;
    rtb_KPCoM_idx_0 =
        torqueBalancingYoga_B.SFunction_j[12] - torqueBalancingYoga_B.SFunction_dk[12];
    rtb_KDCoM_idx_0 =
        torqueBalancingYoga_B.SFunction_oe[12] - torqueBalancingYoga_B.SFunction_dk[12];
    gravityWrench[4] = 0.0;
    rtb_KPCoM_idx_1 =
        torqueBalancingYoga_B.SFunction_j[13] - torqueBalancingYoga_B.SFunction_dk[13];
    rtb_KDCoM_idx_1 =
        torqueBalancingYoga_B.SFunction_oe[13] - torqueBalancingYoga_B.SFunction_dk[13];
    gravityWrench[5] = 0.0;
    rtb_KPCoM_idx_2 =
        torqueBalancingYoga_B.SFunction_j[14] - torqueBalancingYoga_B.SFunction_dk[14];
    rtb_KDCoM_idx_2 =
        torqueBalancingYoga_B.SFunction_oe[14] - torqueBalancingYoga_B.SFunction_dk[14];
    for (i_0 = 0; i_0 < 6; i_0++) {
        pinvJb[6 * i_0] = d[3 * i_0];
        pinvJb[1 + 6 * i_0] = d[3 * i_0 + 1];
        pinvJb[2 + 6 * i_0] = d[3 * i_0 + 2];
    }

    pinvJb[3] = 0.0;
    pinvJb[9] = -rtb_KDCoM_idx_2;
    pinvJb[15] = rtb_KDCoM_idx_1;
    pinvJb[4] = rtb_KDCoM_idx_2;
    pinvJb[10] = 0.0;
    pinvJb[16] = -rtb_KDCoM_idx_0;
    pinvJb[5] = -rtb_KDCoM_idx_1;
    pinvJb[11] = rtb_KDCoM_idx_0;
    pinvJb[17] = 0.0;
    for (i_0 = 0; i_0 < 3; i_0++) {
        Gamma_tmp = 6 * (i_0 + 3);
        pinvJb[3 + Gamma_tmp] = c[3 * i_0];
        pinvJb[4 + Gamma_tmp] = c[3 * i_0 + 1];
        pinvJb[5 + Gamma_tmp] = c[3 * i_0 + 2];
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        AR[6 * i_0] = d[3 * i_0];
        AR[1 + 6 * i_0] = d[3 * i_0 + 1];
        AR[2 + 6 * i_0] = d[3 * i_0 + 2];
    }

    AR[3] = 0.0;
    AR[9] = -rtb_KPCoM_idx_2;
    AR[15] = rtb_KPCoM_idx_1;
    AR[4] = rtb_KPCoM_idx_2;
    AR[10] = 0.0;
    AR[16] = -rtb_KPCoM_idx_0;
    AR[5] = -rtb_KPCoM_idx_1;
    AR[11] = rtb_KPCoM_idx_0;
    AR[17] = 0.0;
    for (i_0 = 0; i_0 < 3; i_0++) {
        AR_tmp = 6 * (i_0 + 3);
        AR[3 + AR_tmp] = c[3 * i_0];
        AR[4 + AR_tmp] = c[3 * i_0 + 1];
        AR[5 + AR_tmp] = c[3 * i_0 + 2];
    }

    for (i_0 = 0; i_0 < 36; i_0++) {
        pinvJb_0[i_0] = pinvJb[i_0];
        pinvJb_0[i_0 + 36] = AR[i_0];
    }

    torqueBalancingYoga_inv(pinvJb, pinvJb_3);
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            i = AR_tmp + 12 * i_0;
            a[i] = pinvJb_3[6 * i_0 + AR_tmp];
            a[i + 6] = 0.0;
        }
    }

    torqueBalancingYoga_inv(AR, pinvJb_3);
    torqueBalancingYoga_pinv(pinvJb_0, torqueBalancingYoga_P.Reg.pinvTol, pinvJb_1);
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            i = AR_tmp + 12 * i_0;
            tmp_2[i] = 0.0;
            tmp_2[i + 6] = pinvJb_3[6 * i_0 + AR_tmp] * torqueBalancingYoga_B.MultiportSwitch1[33]
                           * (1.0 - torqueBalancingYoga_B.MultiportSwitch1[32]);
        }

        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            i = 12 * i_0 + AR_tmp;
            pinvA[AR_tmp + 12 * i_0] = (pinvJb_1[i] * torqueBalancingYoga_B.MultiportSwitch1[32]
                                            * torqueBalancingYoga_B.MultiportSwitch1[33]
                                        + a[i] * torqueBalancingYoga_B.MultiportSwitch1[32]
                                              * (1.0 - torqueBalancingYoga_B.MultiportSwitch1[33]))
                                       + tmp_2[i];
        }
    }

    torqueBalancingYoga_eye_a(JcMinv_1);
    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 6; i++) {
                acc_CoM_des_idx_1 += pinvA[12 * i + AR_tmp] * pinvJb_0[6 * i_0 + i];
            }

            tmp_3[AR_tmp + 12 * i_0] = JcMinv_1[12 * i_0 + AR_tmp] - acc_CoM_des_idx_1;
        }
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            rtb_NA[AR_tmp + 12 * i_0] = tmp_3[12 * i_0 + AR_tmp]
                                        * torqueBalancingYoga_B.MultiportSwitch1[32]
                                        * torqueBalancingYoga_B.MultiportSwitch1[33];
        }
    }

    for (i_0 = 0; i_0 < 29; i_0++) {
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            Gamma_tmp = 6 * i_0 + AR_tmp;
            i = AR_tmp + 12 * i_0;
            Jc[i] = torqueBalancingYoga_B.SFunction_c4[Gamma_tmp]
                    * torqueBalancingYoga_B.MultiportSwitch1[32];
            Jc[i + 6] = torqueBalancingYoga_B.SFunction_b3[Gamma_tmp]
                        * torqueBalancingYoga_B.MultiportSwitch1[33];
        }
    }

    torqueBalancingYoga_mrdivide_h(Jc, torqueBalancingYoga_B.M_with_inertia, JcMinv);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            Gamma_tmp = AR_tmp + 12 * i_0;
            JcMinvSt[Gamma_tmp] = 0.0;
            for (i = 0; i < 29; i++) {
                JcMinvSt[Gamma_tmp] = JcMinv[12 * i + AR_tmp] * (real_T) c_b[29 * i_0 + i]
                                      + JcMinvSt[12 * i_0 + AR_tmp];
            }
        }
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            pinvJb_2[AR_tmp + 23 * i_0] =
                torqueBalancingYoga_B.M_with_inertia[(6 + AR_tmp) * 29 + i_0];
        }

        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            pinvJb_3[AR_tmp + 6 * i_0] = torqueBalancingYoga_B.M_with_inertia[29 * i_0 + AR_tmp];
        }
    }

    torqueBalancingYoga_mrdivide_hg(pinvJb_2, pinvJb_3);
    torqueBalancingYoga_pinvDamped(JcMinvSt, torqueBalancingYoga_P.Reg.pinvDamp, Pinv_JcMinvSt);
    torqueBalancingYoga_eye_a0(torqueBalancingYoga_B.Gamma_m);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 12; i++) {
                acc_CoM_des_idx_1 += Pinv_JcMinvSt[23 * i + AR_tmp] * JcMinvSt[12 * i_0 + i];
            }

            torqueBalancingYoga_B.invTGamma[AR_tmp + 23 * i_0] =
                torqueBalancingYoga_B.Gamma_m[23 * i_0 + AR_tmp] - acc_CoM_des_idx_1;
        }
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            c_a[AR_tmp + 23 * i_0] = torqueBalancingYoga_B.M_with_inertia[(6 + AR_tmp) * 29 + i_0];
        }

        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            pinvJb_3[AR_tmp + 6 * i_0] = torqueBalancingYoga_B.M_with_inertia[29 * i_0 + AR_tmp];
        }
    }

    torqueBalancingYoga_mrdivide_hg(c_a, pinvJb_3);
    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 6; i++) {
                acc_CoM_des_idx_1 +=
                    torqueBalancingYoga_B.M_with_inertia[(6 + AR_tmp) * 29 + i] * c_a[23 * i + i_0];
            }

            torqueBalancingYoga_B.Gamma_m[i_0 + 23 * AR_tmp] =
                torqueBalancingYoga_B.M_with_inertia[((6 + AR_tmp) * 29 + i_0) + 6]
                - acc_CoM_des_idx_1;
        }
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            Gamma_tmp = AR_tmp + 23 * i_0;
            torqueBalancingYoga_B.Gamma[Gamma_tmp] = 0.0;
            for (i = 0; i < 23; i++) {
                torqueBalancingYoga_B.Gamma[Gamma_tmp] =
                    torqueBalancingYoga_B.invTGamma[23 * i + AR_tmp]
                        * torqueBalancingYoga_B.Gamma_m[23 * i_0 + i]
                    + torqueBalancingYoga_B.Gamma[23 * i_0 + AR_tmp];
            }
        }
    }

    for (i_0 = 0; i_0 < 3; i_0++) {
        tmp[3 * i_0] = torqueBalancingYoga_B.SFunction_oe[i_0];
        tmp_0[3 * i_0] = torqueBalancingYoga_B.SFunction_oe[i_0];
        AR_tmp = 1 + 3 * i_0;
        tmp[AR_tmp] = torqueBalancingYoga_B.SFunction_oe[i_0 + 4];
        tmp_0[AR_tmp] = torqueBalancingYoga_B.SFunction_oe[i_0 + 4];
        AR_tmp = 2 + 3 * i_0;
        tmp[AR_tmp] = torqueBalancingYoga_B.SFunction_oe[i_0 + 8];
        tmp_0[AR_tmp] = torqueBalancingYoga_B.SFunction_oe[i_0 + 8];
    }

    torqueBalancingYoga_blkdiag(tmp, tmp_0, pinvJb_3);
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 19; AR_tmp++) {
            Gamma_tmp = AR_tmp + 19 * i_0;
            torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[Gamma_tmp] = 0.0;
            for (i = 0; i < 6; i++) {
                torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[Gamma_tmp] =
                    torqueBalancingYoga_P.ConstraintsMatrix[19 * i + AR_tmp] * pinvJb_3[6 * i_0 + i]
                    + torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[19 * i_0 + AR_tmp];
            }
        }
    }

    for (i_0 = 0; i_0 < 3; i_0++) {
        tmp[3 * i_0] = torqueBalancingYoga_B.SFunction_j[i_0];
        tmp_0[3 * i_0] = torqueBalancingYoga_B.SFunction_j[i_0];
        AR_tmp = 1 + 3 * i_0;
        tmp[AR_tmp] = torqueBalancingYoga_B.SFunction_j[i_0 + 4];
        tmp_0[AR_tmp] = torqueBalancingYoga_B.SFunction_j[i_0 + 4];
        AR_tmp = 2 + 3 * i_0;
        tmp[AR_tmp] = torqueBalancingYoga_B.SFunction_j[i_0 + 8];
        tmp_0[AR_tmp] = torqueBalancingYoga_B.SFunction_j[i_0 + 8];
    }

    torqueBalancingYoga_blkdiag(tmp, tmp_0, pinvJb_3);
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 19; AR_tmp++) {
            Gamma_tmp = AR_tmp + 19 * i_0;
            constraintMatrixRightFoot[Gamma_tmp] = 0.0;
            for (i = 0; i < 6; i++) {
                constraintMatrixRightFoot[Gamma_tmp] =
                    torqueBalancingYoga_P.ConstraintsMatrix[19 * i + AR_tmp] * pinvJb_3[6 * i_0 + i]
                    + constraintMatrixRightFoot[19 * i_0 + AR_tmp];
            }
        }
    }

    torqueBalancingYoga_blkdiag_a(torqueBalancingYoga_B.ConstraintsMatrixQP1Foot,
                                  constraintMatrixRightFoot,
                                  ConstraintsMatrix2Feet);
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            c_a[AR_tmp + 23 * i_0] = torqueBalancingYoga_B.M_with_inertia[(6 + AR_tmp) * 29 + i_0];
        }

        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            pinvJb_3[AR_tmp + 6 * i_0] = torqueBalancingYoga_B.M_with_inertia[29 * i_0 + AR_tmp];
        }
    }

    torqueBalancingYoga_mrdivide_hg(c_a, pinvJb_3);
    torqueBalancingYoga_diag(&torqueBalancingYoga_B.MultiportSwitch1[34],
                             torqueBalancingYoga_B.Gamma_m);
    torqueBalancingYoga_pinv_o(
        torqueBalancingYoga_B.Gamma, torqueBalancingYoga_P.Reg.pinvTol, torqueBalancingYoga_B.dv0);
    torqueBalancingYoga_diag(torqueBalancingYoga_P.Gain.dampings, torqueBalancingYoga_B.dv1);
    torqueBalancingYoga_pinv_o(
        torqueBalancingYoga_B.Gamma, torqueBalancingYoga_P.Reg.pinvTol, torqueBalancingYoga_B.dv2);
    for (i_0 = 0; i_0 < 12; i_0++) {
        rtb_fLdotDesC1C2[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 29; AR_tmp++) {
            rtb_fLdotDesC1C2[i_0] +=
                JcMinv[12 * AR_tmp + i_0] * torqueBalancingYoga_B.SFunction_f[AR_tmp];
        }
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        pinvA_0[i_0] =
            torqueBalancingYoga_B.SFunction_ej[i_0] * torqueBalancingYoga_B.MultiportSwitch1[32];
        pinvA_0[i_0 + 6] =
            torqueBalancingYoga_B.SFunction_h[i_0] * torqueBalancingYoga_B.MultiportSwitch1[33];
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        JcMinv_0[i_0] = rtb_fLdotDesC1C2[i_0] - pinvA_0[i_0];
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 23; i++) {
                acc_CoM_des_idx_1 += torqueBalancingYoga_B.Gamma_m[23 * i + i_0]
                                     * torqueBalancingYoga_B.dv0[23 * AR_tmp + i];
            }

            torqueBalancingYoga_B.dv3[i_0 + 23 * AR_tmp] =
                (real_T) d_b[23 * AR_tmp + i_0] * torqueBalancingYoga_P.Reg.impedances
                + acc_CoM_des_idx_1;
        }

        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            i = i_0 + 23 * AR_tmp;
            torqueBalancingYoga_B.dv4[i] = 0.0;
            for (Gamma_tmp = 0; Gamma_tmp < 23; Gamma_tmp++) {
                torqueBalancingYoga_B.dv4[i] =
                    torqueBalancingYoga_B.dv3[23 * Gamma_tmp + i_0]
                        * torqueBalancingYoga_B.Gamma[23 * AR_tmp + Gamma_tmp]
                    + torqueBalancingYoga_B.dv4[23 * AR_tmp + i_0];
            }
        }

        invTGamma[i_0] =
            torqueBalancingYoga_B.MultiportSwitch1[87 + i_0] - torqueBalancingYoga_B.Switch5[i_0];
        acc_CoM_des_idx_1 = 0.0;
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            acc_CoM_des_idx_1 += c_a[23 * AR_tmp + i_0] * torqueBalancingYoga_B.SFunction_f[AR_tmp];
        }

        tmp_4[i_0] = torqueBalancingYoga_B.SFunction_f[6 + i_0] - acc_CoM_des_idx_1;
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        rtb_impedances[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 23; i++) {
                acc_CoM_des_idx_1 += torqueBalancingYoga_B.dv1[23 * i + i_0]
                                     * torqueBalancingYoga_B.dv2[23 * AR_tmp + i];
            }

            i = 23 * AR_tmp + i_0;
            torqueBalancingYoga_B.Gamma_m[i_0 + 23 * AR_tmp] =
                (real_T) d_b[i] * torqueBalancingYoga_P.Reg.dampings + acc_CoM_des_idx_1;
            rtb_impedances[i_0] += torqueBalancingYoga_B.dv4[i] * invTGamma[AR_tmp];
        }

        acc_CoM_des_idx_1 = 0.0;
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            i = i_0 + 23 * AR_tmp;
            torqueBalancingYoga_B.dv0[i] = 0.0;
            for (Gamma_tmp = 0; Gamma_tmp < 23; Gamma_tmp++) {
                torqueBalancingYoga_B.dv0[i] =
                    torqueBalancingYoga_B.Gamma_m[23 * Gamma_tmp + i_0]
                        * torqueBalancingYoga_B.Gamma[23 * AR_tmp + Gamma_tmp]
                    + torqueBalancingYoga_B.dv0[23 * AR_tmp + i_0];
            }

            acc_CoM_des_idx_1 += torqueBalancingYoga_B.dv0[23 * AR_tmp + i_0]
                                 * torqueBalancingYoga_B.MultiportSwitch1[64 + AR_tmp];
        }

        rtb_impedances_0[i_0] = (tmp_4[i_0] - rtb_impedances[i_0]) - acc_CoM_des_idx_1;
        Pinv_JcMinvSt_0[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            Pinv_JcMinvSt_0[i_0] += Pinv_JcMinvSt[23 * AR_tmp + i_0] * JcMinv_0[AR_tmp];
        }
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        invTGamma[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            invTGamma[i_0] +=
                torqueBalancingYoga_B.invTGamma[23 * AR_tmp + i_0] * rtb_impedances_0[AR_tmp];
        }

        rtb_impedances[i_0] = Pinv_JcMinvSt_0[i_0] + invTGamma[i_0];
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            Gamma_tmp = i_0 + 12 * AR_tmp;
            JcMinv_1[Gamma_tmp] = 0.0;
            for (i = 0; i < 29; i++) {
                JcMinv_1[Gamma_tmp] =
                    JcMinv[12 * i + i_0] * Jc[12 * i + AR_tmp] + JcMinv_1[12 * AR_tmp + i_0];
            }
        }
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 6; i++) {
                acc_CoM_des_idx_1 += pinvJb_2[23 * i + i_0] * Jc[12 * i + AR_tmp];
            }

            Gamma_tmp = i_0 + 23 * AR_tmp;
            rtb_Sigma[Gamma_tmp] = Jc[(6 + i_0) * 12 + AR_tmp] - acc_CoM_des_idx_1;
            Pinv_JcMinvSt_1[Gamma_tmp] = 0.0;
            for (i = 0; i < 12; i++) {
                Pinv_JcMinvSt_1[Gamma_tmp] = Pinv_JcMinvSt[23 * i + i_0] * JcMinv_1[12 * AR_tmp + i]
                                             + Pinv_JcMinvSt_1[23 * AR_tmp + i_0];
            }
        }
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            Gamma_tmp = AR_tmp + 23 * i_0;
            Pinv_JcMinvSt[Gamma_tmp] = 0.0;
            for (i = 0; i < 23; i++) {
                Pinv_JcMinvSt[Gamma_tmp] =
                    torqueBalancingYoga_B.invTGamma[23 * i + AR_tmp] * rtb_Sigma[23 * i_0 + i]
                    + Pinv_JcMinvSt[23 * i_0 + AR_tmp];
            }
        }
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            i = 23 * i_0 + AR_tmp;
            rtb_Sigma[AR_tmp + 23 * i_0] = -(Pinv_JcMinvSt_1[i] + Pinv_JcMinvSt[i]);
        }
    }

    for (i_0 = 0; i_0 < 3; i_0++) {
        acc_CoM_des_idx_1 = 0.0;
        for (AR_tmp = 0; AR_tmp < 29; AR_tmp++) {
            acc_CoM_des_idx_1 += torqueBalancingYoga_B.SFunction_n[6 * AR_tmp + i_0]
                                 * torqueBalancingYoga_B.MultiportSwitch1[58 + AR_tmp];
        }

        LDotDes[i_0] = ((rtb_Switch3[6 + i_0]
                         - (torqueBalancingYoga_B.SFunction_dk[12 + i_0] - rtb_Switch3[i_0])
                               * torqueBalancingYoga_B.MultiportSwitch1[111 + i_0])
                        - (acc_CoM_des_idx_1 - rtb_Switch3[3 + i_0])
                              * torqueBalancingYoga_B.MultiportSwitch1[114 + i_0])
                       * torqueBalancingYoga_B.M_with_inertia[0];
        LDotDes[i_0 + 3] = torqueBalancingYoga_B.SFunction_c[3 + i_0]
                               * -torqueBalancingYoga_P.Gain.KD_AngularMomentum
                           - torqueBalancingYoga_B.SFunction_on[3 + i_0]
                                 * torqueBalancingYoga_P.Gain.KP_AngularMomentum;
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        LDotDes_0[i_0] = LDotDes[i_0] - gravityWrench[i_0];
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        pinvA_0[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            pinvA_0[i_0] += pinvA[12 * AR_tmp + i_0] * LDotDes_0[AR_tmp];
        }

        rtb_fLdotDesC1C2[i_0] = pinvA_0[i_0] * torqueBalancingYoga_B.MultiportSwitch1[32]
                                * torqueBalancingYoga_B.MultiportSwitch1[33];
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            Gamma_tmp = AR_tmp + 23 * i_0;
            JcMinvSt[Gamma_tmp] = 0.0;
            for (i = 0; i < 12; i++) {
                JcMinvSt[Gamma_tmp] =
                    rtb_Sigma[23 * i + AR_tmp] * rtb_NA[12 * i_0 + i] + JcMinvSt[23 * i_0 + AR_tmp];
            }
        }
    }

    rtb_Clock1 = (1.0 - torqueBalancingYoga_B.MultiportSwitch1[33])
                 * torqueBalancingYoga_B.MultiportSwitch1[32];
    acc_CoM_des_idx_0 = (1.0 - torqueBalancingYoga_B.MultiportSwitch1[32])
                        * torqueBalancingYoga_B.MultiportSwitch1[33];
    for (i_0 = 0; i_0 < 36; i_0++) {
        pinvJb[i_0] = pinvJb[i_0] * torqueBalancingYoga_B.MultiportSwitch1[32]
                          * (1.0 - torqueBalancingYoga_B.MultiportSwitch1[33])
                      + AR[i_0] * torqueBalancingYoga_B.MultiportSwitch1[33]
                            * (1.0 - torqueBalancingYoga_B.MultiportSwitch1[32]);
    }

    torqueBalancingYog_pinvDamped_a(
        JcMinvSt, torqueBalancingYoga_P.Reg.pinvDamp * 1.0E-5, unusedExpr);

    /* '<S108>:1:9' */
    /* '<S108>:1:10' */
    /* '<S108>:1:11' */
    /* '<S108>:1:12' */
    for (i_0 = 0; i_0 < 6; i_0++) {
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 6; i++) {
                acc_CoM_des_idx_1 += pinvJb[6 * i_0 + i] * pinvJb[6 * AR_tmp + i];
            }

            i = 6 * AR_tmp + i_0;
            torqueBalancingYoga_B.HessianMatrixQP1Foot[i_0 + 6 * AR_tmp] =
                (real_T) f_a[i] * torqueBalancingYoga_P.Reg.HessianQP + acc_CoM_des_idx_1;
            pinvJb_3[AR_tmp + 6 * i_0] = -pinvJb[i];
        }

        LDotDes_0[i_0] = LDotDes[i_0] - gravityWrench[i_0];
    }

    for (i_0 = 0; i_0 < 6; i_0++) {
        torqueBalancingYoga_B.gradientQP1Foot[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
            torqueBalancingYoga_B.gradientQP1Foot[i_0] +=
                pinvJb_3[6 * AR_tmp + i_0] * LDotDes_0[AR_tmp];
        }
    }

    for (i_0 = 0; i_0 < 114; i_0++) {
        torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[i_0] =
            rtb_Clock1 * torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[i_0]
            + acc_CoM_des_idx_0 * constraintMatrixRightFoot[i_0];
    }

    memcpy(&torqueBalancingYoga_B.bVectorConstraintsQp1Foot[0],
           &torqueBalancingYoga_P.bVectorConstraints[0],
           19U * sizeof(real_T));
    for (i_0 = 0; i_0 < 12; i_0++) {
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            acc_CoM_des_idx_1 = 0.0;
            for (i = 0; i < 23; i++) {
                acc_CoM_des_idx_1 += JcMinvSt[23 * i_0 + i] * JcMinvSt[23 * AR_tmp + i];
            }

            torqueBalancingYoga_B.HessianMatrixQP2Feet[i_0 + 12 * AR_tmp] =
                (real_T) g_a[12 * AR_tmp + i_0] * torqueBalancingYoga_P.Reg.HessianQP
                + acc_CoM_des_idx_1;
        }
    }

    for (i_0 = 0; i_0 < 23; i_0++) {
        acc_CoM_des_idx_1 = 0.0;
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            acc_CoM_des_idx_1 += rtb_Sigma[23 * AR_tmp + i_0] * rtb_fLdotDesC1C2[AR_tmp];
        }

        rtb_impedances_0[i_0] = rtb_impedances[i_0] + acc_CoM_des_idx_1;
    }

    for (i_0 = 0; i_0 < 12; i_0++) {
        torqueBalancingYoga_B.gradientQP2Feet[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
            torqueBalancingYoga_B.gradientQP2Feet[i_0] +=
                JcMinvSt[23 * i_0 + AR_tmp] * rtb_impedances_0[AR_tmp];
        }

        for (AR_tmp = 0; AR_tmp < 38; AR_tmp++) {
            Gamma_tmp = AR_tmp + 38 * i_0;
            torqueBalancingYoga_B.ConstraintsMatrixQP2Feet[Gamma_tmp] = 0.0;
            for (i = 0; i < 12; i++) {
                torqueBalancingYoga_B.ConstraintsMatrixQP2Feet[Gamma_tmp] =
                    ConstraintsMatrix2Feet[38 * i + AR_tmp] * rtb_NA[12 * i_0 + i]
                    + torqueBalancingYoga_B.ConstraintsMatrixQP2Feet[38 * i_0 + AR_tmp];
            }
        }
    }

    for (i_0 = 0; i_0 < 19; i_0++) {
        tmp_5[i_0] = torqueBalancingYoga_P.bVectorConstraints[i_0];
        tmp_5[i_0 + 19] = torqueBalancingYoga_P.bVectorConstraints[i_0];
    }

    for (i_0 = 0; i_0 < 38; i_0++) {
        ConstraintsMatrix2Feet_0[i_0] = 0.0;
        for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
            ConstraintsMatrix2Feet_0[i_0] +=
                ConstraintsMatrix2Feet[38 * AR_tmp + i_0] * rtb_fLdotDesC1C2[AR_tmp];
        }

        torqueBalancingYoga_B.bVectorConstraintsQp2Feet[i_0] =
            tmp_5[i_0] - ConstraintsMatrix2Feet_0[i_0];
    }

    /* End of MATLAB Function: '<S6>/Balancing Controller ' */

    /* MATLAB Function: '<S142>/ContactsTransition' */
    /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/ContactsTransition': '<S143>:1'
     */
    /* '<S143>:1:3' */
    rtb_Clock1 =
        torqueBalancingYoga_B.MultiportSwitch1[32] + torqueBalancingYoga_B.MultiportSwitch1[33];
    torqueBalancingYoga_B.onOneFoot = ((!(rtb_Clock1 > 1.9)) && (rtb_Clock1 > 0.9));

    /* SignalConversion: '<S142>/HiddenBuf_InsertedFor_One Foot_at_inport_4' */
    torqueBalancingYoga_B.HiddenBuf_InsertedFor_OneFoot_a = torqueBalancingYoga_B.onOneFoot;

    /* Outputs for Enabled SubSystem: '<S142>/One Foot' incorporates:
     *  EnablePort: '<S144>/Enable'
     */
    if (torqueBalancingYoga_B.HiddenBuf_InsertedFor_OneFoot_a) {
        if (!torqueBalancingYoga_DW.OneFoot_MODE) {
            torqueBalancingYoga_DW.OneFoot_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.OneFoot_MODE) {
            torqueBalancingYoga_DW.OneFoot_MODE = false;
        }
    }

    if (torqueBalancingYoga_DW.OneFoot_MODE) {
        /* S-Function (WBToolbox): '<S144>/QP Two Feet' */
        /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/One Foot/Analytical
         * Solution One Foot (unconstrained)': '<S148>:1' */
        /* '<S148>:1:3' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S144>/QP Two Feet

        /* MATLAB Function: '<S144>/Process QP output' incorporates:
         *  MATLAB Function: '<S144>/Analytical Solution One Foot (unconstrained)'
         */
        /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/One Foot/Process QP
         * output': '<S150>:1' */
        /* '<S150>:1:3' */
        if (torqueBalancingYoga_P.Config.USE_QP_SOLVER
            && (std::abs(torqueBalancingYoga_B.QPTwoFeet_o2_h) < 0.01)) {
            for (i = 0; i < 6; i++) {
                torqueBalancingYoga_B.f0_h[i] = torqueBalancingYoga_B.QPTwoFeet_o1_d[i];
            }
        }
        else {
            /* MATLAB Function: '<S144>/Analytical Solution One Foot (unconstrained)' */
            torqueBalancingYoga_invNxN(torqueBalancingYoga_B.HessianMatrixQP1Foot, pinvJb_3);
            for (i_0 = 0; i_0 < 6; i_0++) {
                for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                    AR[AR_tmp + 6 * i_0] = -pinvJb_3[6 * i_0 + AR_tmp];
                }
            }

            for (i_0 = 0; i_0 < 6; i_0++) {
                torqueBalancingYoga_B.f0_h[i_0] = 0.0;
                for (AR_tmp = 0; AR_tmp < 6; AR_tmp++) {
                    torqueBalancingYoga_B.f0_h[i_0] +=
                        AR[6 * AR_tmp + i_0] * torqueBalancingYoga_B.gradientQP1Foot[AR_tmp];
                }
            }
        }

        /* End of MATLAB Function: '<S144>/Process QP output' */
    }

    /* End of Outputs for SubSystem: '<S142>/One Foot' */

    /* MATLAB Function: '<S142>/Process One Foot Output' */
    /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/Process One Foot Output':
     * '<S145>:1' */
    /* '<S145>:1:3' */
    rtb_Clock1 = std::abs(torqueBalancingYoga_B.MultiportSwitch1[32]
                          - torqueBalancingYoga_B.MultiportSwitch1[33]);

    /* Logic: '<S142>/not' */
    torqueBalancingYoga_B.not_n = !torqueBalancingYoga_B.onOneFoot;

    /* SignalConversion: '<S142>/HiddenBuf_InsertedFor_Two Feet_at_inport_4' */
    torqueBalancingYoga_B.HiddenBuf_InsertedFor_TwoFeet_a = torqueBalancingYoga_B.not_n;

    /* Outputs for Enabled SubSystem: '<S142>/Two Feet' incorporates:
     *  EnablePort: '<S146>/Enable'
     */
    if (torqueBalancingYoga_B.HiddenBuf_InsertedFor_TwoFeet_a) {
        if (!torqueBalancingYoga_DW.TwoFeet_MODE) {
            torqueBalancingYoga_DW.TwoFeet_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.TwoFeet_MODE) {
            torqueBalancingYoga_DW.TwoFeet_MODE = false;
        }
    }

    if (torqueBalancingYoga_DW.TwoFeet_MODE) {
        /* S-Function (WBToolbox): '<S146>/QP Two Feet' */
        /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/Two Feet/Analytical
         * Solution Two Feet (unconstrained)': '<S151>:1' */
        /* '<S151>:1:3' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S146>/QP Two Feet

        /* MATLAB Function: '<S146>/Process QP output' incorporates:
         *  MATLAB Function: '<S146>/Analytical Solution Two Feet (unconstrained)'
         */
        /* MATLAB Function 'controller_QP/Compute joint torques/QPSolver/Two Feet/Process QP
         * output': '<S153>:1' */
        /* '<S153>:1:3' */
        if (torqueBalancingYoga_P.Config.USE_QP_SOLVER
            && (std::abs(torqueBalancingYoga_B.QPTwoFeet_o2) < 0.01)) {
            memcpy(&torqueBalancingYoga_B.f0[0],
                   &torqueBalancingYoga_B.QPTwoFeet_o1[0],
                   12U * sizeof(real_T));
        }
        else {
            /* MATLAB Function: '<S146>/Analytical Solution Two Feet (unconstrained)' */
            torqueBalancingYoga_invNxN_l(torqueBalancingYoga_B.HessianMatrixQP2Feet, JcMinv_1);
            for (i_0 = 0; i_0 < 12; i_0++) {
                for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                    tmp_3[AR_tmp + 12 * i_0] = -JcMinv_1[12 * i_0 + AR_tmp];
                }
            }

            for (i_0 = 0; i_0 < 12; i_0++) {
                torqueBalancingYoga_B.f0[i_0] = 0.0;
                for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                    torqueBalancingYoga_B.f0[i_0] +=
                        tmp_3[12 * AR_tmp + i_0] * torqueBalancingYoga_B.gradientQP2Feet[AR_tmp];
                }
            }
        }

        /* End of MATLAB Function: '<S146>/Process QP output' */
    }

    /* End of Outputs for SubSystem: '<S142>/Two Feet' */

    /* Switch: '<S111>/Switch' incorporates:
     *  Constant: '<S111>/ '
     *  Product: '<S110>/Product1'
     *  Sum: '<S110>/Add'
     */
    if (torqueBalancingYoga_P.Config.USE_MOTOR_REFLECTED_INERTIA) {
        /* MATLAB Function: '<S142>/Process One Foot Output' incorporates:
         *  Sum: '<S110>/Sum'
         */
        for (i_0 = 0; i_0 < 6; i_0++) {
            pinvA_0[i_0] = torqueBalancingYoga_B.f0_h[i_0]
                           * torqueBalancingYoga_B.MultiportSwitch1[32] * rtb_Clock1;
            pinvA_0[i_0 + 6] = torqueBalancingYoga_B.f0_h[i_0]
                               * torqueBalancingYoga_B.MultiportSwitch1[33] * rtb_Clock1;
        }

        /* Sum: '<S110>/Add1' incorporates:
         *  Product: '<S110>/Product1'
         *  Product: '<S110>/Product2'
         *  Sum: '<S110>/Sum'
         */
        for (i_0 = 0; i_0 < 12; i_0++) {
            acc_CoM_des_idx_1 = 0.0;
            for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                acc_CoM_des_idx_1 += rtb_NA[12 * AR_tmp + i_0] * torqueBalancingYoga_B.f0[AR_tmp];
            }

            JcMinv_0[i_0] = (rtb_fLdotDesC1C2[i_0] + pinvA_0[i_0]) + acc_CoM_des_idx_1;
        }

        /* Product: '<S111>/Product' incorporates:
         *  Gain: '<S111>/Gain'
         *  Product: '<S110>/Product1'
         *  Sum: '<S110>/Add'
         *  Sum: '<S111>/Add'
         */
        for (i_0 = 0; i_0 < 23; i_0++) {
            invTGamma[i_0] = 0.0;
            for (AR_tmp = 0; AR_tmp < 23; AR_tmp++) {
                invTGamma[i_0] += torqueBalancingYoga_B.reflectedInertia[23 * AR_tmp + i_0]
                                  * torqueBalancingYoga_B.SFunction[AR_tmp];
            }

            acc_CoM_des_idx_1 = 0.0;
            for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                acc_CoM_des_idx_1 += rtb_Sigma[23 * AR_tmp + i_0] * JcMinv_0[AR_tmp];
            }

            rtb_Switch_c[i_0] = (rtb_impedances[i_0] + acc_CoM_des_idx_1)
                                - torqueBalancingYoga_P.Config.K_ff * invTGamma[i_0];
        }

        /* End of Product: '<S111>/Product' */
    }
    else {
        /* MATLAB Function: '<S142>/Process One Foot Output' incorporates:
         *  Sum: '<S110>/Sum'
         */
        for (i_0 = 0; i_0 < 6; i_0++) {
            pinvA_0[i_0] = torqueBalancingYoga_B.f0_h[i_0]
                           * torqueBalancingYoga_B.MultiportSwitch1[32] * rtb_Clock1;
            pinvA_0[i_0 + 6] = torqueBalancingYoga_B.f0_h[i_0]
                               * torqueBalancingYoga_B.MultiportSwitch1[33] * rtb_Clock1;
        }

        /* Sum: '<S110>/Add1' incorporates:
         *  Product: '<S110>/Product1'
         *  Product: '<S110>/Product2'
         *  Sum: '<S110>/Sum'
         */
        for (i_0 = 0; i_0 < 12; i_0++) {
            acc_CoM_des_idx_1 = 0.0;
            for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                acc_CoM_des_idx_1 += rtb_NA[12 * AR_tmp + i_0] * torqueBalancingYoga_B.f0[AR_tmp];
            }

            JcMinv_0[i_0] = (rtb_fLdotDesC1C2[i_0] + pinvA_0[i_0]) + acc_CoM_des_idx_1;
        }

        for (i_0 = 0; i_0 < 23; i_0++) {
            acc_CoM_des_idx_1 = 0.0;
            for (AR_tmp = 0; AR_tmp < 12; AR_tmp++) {
                acc_CoM_des_idx_1 += rtb_Sigma[23 * AR_tmp + i_0] * JcMinv_0[AR_tmp];
            }

            rtb_Switch_c[i_0] = rtb_impedances[i_0] + acc_CoM_des_idx_1;
        }
    }

    /* End of Switch: '<S111>/Switch' */

    /* RelationalOperator: '<S163>/Compare' incorporates:
     *  Clock: '<S162>/Clock'
     *  Constant: '<S163>/Constant'
     */
    rtb_Compare_a =
        ((&torqueBalancingYoga_M)->Timing.t[0] == torqueBalancingYoga_P.CompareToConstant_const_l);

    /* MATLAB Function: '<S162>/MATLAB Function' */
    torqueBalancin_MATLABFunction_k(rtb_Switch_c,
                                    &torqueBalancingYoga_B.sf_MATLABFunction_p,
                                    &torqueBalancingYoga_DW.sf_MATLABFunction_p);

    /* MATLAB Function: '<S9>/Saturate the Torque Derivative' */
    /* MATLAB Function 'tauDot Saturation/Saturate the Torque Derivative': '<S161>:1' */
    /* '<S161>:1:3' */
    if (!torqueBalancingYoga_DW.uPrev_not_empty) {
        memcpy(&torqueBalancingYoga_DW.uPrev[0],
               &torqueBalancingYoga_B.sf_MATLABFunction_p.s0[0],
               23U * sizeof(real_T));
        torqueBalancingYoga_DW.uPrev_not_empty = true;
    }

    rtb_Clock1 = torqueBalancingYoga_P.Config.tauDot_maxAbs * torqueBalancingYoga_P.Config.Ts;
    acc_CoM_des_idx_0 =
        -torqueBalancingYoga_P.Config.tauDot_maxAbs * torqueBalancingYoga_P.Config.Ts;
    for (i = 0; i < 23; i++) {
        rtb_KPCoM_idx_0 = rtb_Switch_c[i] - torqueBalancingYoga_DW.uPrev[i];
        if (rtb_KPCoM_idx_0 > rtb_Clock1) {
            rtb_KPCoM_idx_0 = rtb_Clock1;
        }

        if (rtb_KPCoM_idx_0 < acc_CoM_des_idx_0) {
            rtb_KPCoM_idx_0 = acc_CoM_des_idx_0;
        }

        rtb_KPCoM_idx_0 += torqueBalancingYoga_DW.uPrev[i];
        torqueBalancingYoga_DW.uPrev[i] = rtb_KPCoM_idx_0;

        /* Switch: '<S9>/Switch' incorporates:
         *  Constant: '<S9>/Constant'
         */
        if (!torqueBalancingYoga_P.Config.SATURATE_TORQUE_DERIVATIVE) {
            rtb_KPCoM_idx_0 = rtb_Switch_c[i];
        }

        /* End of Switch: '<S9>/Switch' */

        /* Saturate: '<Root>/Saturation' */
        if (rtb_KPCoM_idx_0 > torqueBalancingYoga_P.Sat.torque) {
            torqueBalancingYoga_B.Saturation[i] = torqueBalancingYoga_P.Sat.torque;
        }
        else if (rtb_KPCoM_idx_0 < -torqueBalancingYoga_P.Sat.torque) {
            torqueBalancingYoga_B.Saturation[i] = -torqueBalancingYoga_P.Sat.torque;
        }
        else {
            torqueBalancingYoga_B.Saturation[i] = rtb_KPCoM_idx_0;
        }

        /* End of Saturate: '<Root>/Saturation' */
    }

    /* End of MATLAB Function: '<S9>/Saturate the Torque Derivative' */

    /* S-Function (WBToolbox): '<S5>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S5>/S-Function

    /* Logic: '<S33>/Logical Operator1' incorporates:
     *  Constant: '<S33>/ON_GAZEBO 3'
     *  Constant: '<S33>/ON_GAZEBO 4'
     */
    torqueBalancingYoga_B.LogicalOperator1 =
        (torqueBalancingYoga_P.Config.SCOPES_ALL
         || torqueBalancingYoga_P.Config.SCOPES_GAIN_SCHE_INFO);

    /* Outputs for Enabled SubSystem: '<Root>/emergency stop: joint limits' incorporates:
     *  EnablePort: '<S7>/Enable'
     */
    /* Constant: '<Root>/ON_GAZEBO 1' */
    if (torqueBalancingYoga_P.Config.CHECK_LIMITS) {
        if (!torqueBalancingYoga_DW.emergencystopjointlimits_MODE) {
            torqueBalancingYoga_DW.emergencystopjointlimits_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.emergencystopjointlimits_MODE) {
            torqueBalancingYoga_DW.emergencystopjointlimits_MODE = false;
        }
    }

    /* End of Constant: '<Root>/ON_GAZEBO 1' */
    if (torqueBalancingYoga_DW.emergencystopjointlimits_MODE) {
        /* S-Function (WBToolbox): '<S157>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S157>/S-Function

        /* MATLAB Function: '<S7>/MATLAB Function' incorporates:
         *  Constant: '<S7>/index1'
         */
        /* MATLAB Function 'emergency stop: joint limits/MATLAB Function': '<S158>:1' */
        /* '<S158>:1:3' */
        for (i = 0; i < 23; i++) {
            res[i] =
                ((torqueBalancingYoga_B.MultiportSwitch1[i + 87]
                  < torqueBalancingYoga_B.SFunction_o1[i] + torqueBalancingYoga_P.index1_Value)
                 || (torqueBalancingYoga_B.MultiportSwitch1[i + 87]
                     > torqueBalancingYoga_B.SFunction_o2[i] - torqueBalancingYoga_P.index1_Value));
        }

        rtb_Clock1 = res[0];
        for (i = 0; i < 22; i++) {
            rtb_Clock1 += (real_T) res[i + 1];
        }

        torqueBalancingYoga_B.inRange = (rtb_Clock1 == 0.0);

        /* End of MATLAB Function: '<S7>/MATLAB Function' */

        /* Assertion: '<S7>/Assertion' */
        utAssert(torqueBalancingYoga_B.inRange != 0.0);
    }

    /* End of Outputs for SubSystem: '<Root>/emergency stop: joint limits' */

    /* Logic: '<S2>/Logical Operator1' incorporates:
     *  Constant: '<S2>/ON_GAZEBO 3'
     *  Constant: '<S2>/ON_GAZEBO 4'
     */
    torqueBalancingYoga_B.LogicalOperator1_n =
        (torqueBalancingYoga_P.Config.SCOPES_ALL || torqueBalancingYoga_P.Config.SCOPES_MAIN);

    /* Outputs for Enabled SubSystem: '<S2>/Visualizer' incorporates:
     *  EnablePort: '<S10>/Enable'
     */
    if (torqueBalancingYoga_B.LogicalOperator1_n) {
        if (!torqueBalancingYoga_DW.Visualizer_MODE) {
            torqueBalancingYoga_DW.Visualizer_MODE = true;
        }
    }
    else {
        if (torqueBalancingYoga_DW.Visualizer_MODE) {
            torqueBalancingYoga_DW.Visualizer_MODE = false;
        }
    }

    if (torqueBalancingYoga_DW.Visualizer_MODE) {
        /* S-Function (WBToolbox): '<S12>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S12>/S-Function

        /* S-Function (WBToolbox): '<S11>/S-Function' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr =
                static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S11>/S-Function
    }

    /* End of Outputs for SubSystem: '<S2>/Visualizer' */

    /* Outputs for Enabled SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' incorporates:
     *  EnablePort: '<S159>/Enable'
     */
    /* Constant: '<S8>/ON_GAZEBO ' */
    if (torqueBalancingYoga_P.Config.ON_GAZEBO) {
        /* S-Function (WBToolbox): '<S159>/Simulator Synchronizer' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S159>/Simulator Synchronizer
    }

    /* End of Outputs for SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' */

    /* Outputs for Enabled SubSystem: '<S8>/REAL_TIME_SYNC' incorporates:
     *  EnablePort: '<S160>/Enable'
     */
    /* Logic: '<S8>/Logical Operator' incorporates:
     *  Constant: '<S8>/ON_GAZEBO '
     */
    if (!torqueBalancingYoga_P.Config.ON_GAZEBO) {
        /* S-Function (WBToolbox): '<S160>/Real Time Synchronizer' */
        {
            // Get the CoderBlockInformation from the PWork
            wbt::CoderBlockInformation* blockInfo = nullptr;
            blockInfo = static_cast<wbt::CoderBlockInformation*>(
                torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[1]);

            // Get the Block from the PWork
            wbt::Block* blockPtr = nullptr;
            blockPtr = static_cast<wbt::Block*>(
                torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[0]);

            // Calculate the output
            // --------------------
            bool ok;
            ok = blockPtr->output(blockInfo);

            // Report errors
            if (!ok) {
                std::string error = wbt::Log::getSingleton().getErrors();
                error = "[Output]" + error;

                // Trim the message if needed
                if (error.length() >= 1024) {
                    error = error.substr(0, 1024 - 1);
                }

                // This shouldn't happen
                if (getRTM()->errorStatus) {
                    delete getRTM()->errorStatus;
                    getRTM()->errorStatus = nullptr;
                }

                getRTM()->errorStatus = new char[1024];
                sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
                return;
            }
        }

        // End of S-Function Block: <S160>/Real Time Synchronizer
    }

    /* End of Logic: '<S8>/Logical Operator' */
    /* End of Outputs for SubSystem: '<S8>/REAL_TIME_SYNC' */

    /* S-Function (WBToolbox): '<S8>/Yarp Clock' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[0]);

        // Calculate the output
        // --------------------
        bool ok;
        ok = blockPtr->output(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Output]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S8>/Yarp Clock

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&torqueBalancingYoga_M)->Timing.clockTick0)) {
        ++(&torqueBalancingYoga_M)->Timing.clockTickH0;
    }

    (&torqueBalancingYoga_M)->Timing.t[0] =
        (&torqueBalancingYoga_M)->Timing.clockTick0 * (&torqueBalancingYoga_M)->Timing.stepSize0
        + (&torqueBalancingYoga_M)->Timing.clockTickH0 * (&torqueBalancingYoga_M)->Timing.stepSize0
              * 4294967296.0;

    {
        /* Update absolute timer for sample time: [0.01s, 0.0s] */
        /* The "clockTick1" counts the number of times the code of this task has
         * been executed. The resolution of this integer timer is 0.01, which is the step size
         * of the task. Size of "clockTick1" ensures timer will not overflow during the
         * application lifespan selected.
         * Timer of this task consists of two 32 bit unsigned integers.
         * The two integers represent the low bits Timing.clockTick1 and the high bits
         * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
         */
        (&torqueBalancingYoga_M)->Timing.clockTick1++;
        if (!(&torqueBalancingYoga_M)->Timing.clockTick1) {
            (&torqueBalancingYoga_M)->Timing.clockTickH1++;
        }
    }
}

/* Model initialize function */
void torqueBalancingYogaModelClass::initialize()
{
    /* Registration code */

    /* initialize real-time model */
    (void) memset((void*) (&torqueBalancingYoga_M), 0, sizeof(RT_MODEL_torqueBalancingYoga_T));

    {
        /* Setup solver object */
        rtsiSetSimTimeStepPtr(&(&torqueBalancingYoga_M)->solverInfo,
                              &(&torqueBalancingYoga_M)->Timing.simTimeStep);
        rtsiSetTPtr(&(&torqueBalancingYoga_M)->solverInfo, &rtmGetTPtr((&torqueBalancingYoga_M)));
        rtsiSetStepSizePtr(&(&torqueBalancingYoga_M)->solverInfo,
                           &(&torqueBalancingYoga_M)->Timing.stepSize0);
        rtsiSetErrorStatusPtr(&(&torqueBalancingYoga_M)->solverInfo,
                              (&rtmGetErrorStatus((&torqueBalancingYoga_M))));
        rtsiSetRTModelPtr(&(&torqueBalancingYoga_M)->solverInfo, (&torqueBalancingYoga_M));
    }

    rtsiSetSimTimeStep(&(&torqueBalancingYoga_M)->solverInfo, MAJOR_TIME_STEP);
    rtsiSetSolverName(&(&torqueBalancingYoga_M)->solverInfo, "FixedStepDiscrete");
    rtmSetTPtr(getRTM(), &(&torqueBalancingYoga_M)->Timing.tArray[0]);
    (&torqueBalancingYoga_M)->Timing.stepSize0 = 0.01;

    /* block I/O */
    (void) memset(((void*) &torqueBalancingYoga_B), 0, sizeof(B_torqueBalancingYoga_T));

    /* states (dwork) */
    (void) memset((void*) &torqueBalancingYoga_DW, 0, sizeof(DW_torqueBalancingYoga_T));

    /* Start for S-Function (WBToolbox): '<S156>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Acceleration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S156>/S-Function

    /* Start for S-Function (WBToolbox): '<S31>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Position",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_d[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S31>/S-Function

    /* Start for Enabled SubSystem: '<S4>/State Machine Yoga' */

    /* Start for S-Function (WBToolbox): '<S34>/right_foot_wrench' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/wholeBodyDynamics/right_foot/cartesianEndEffectorWrench:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.wR_WBDT[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S34>/right_foot_wrench

    /* Start for S-Function (WBToolbox): '<S34>/left_foot_wrench' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/wholeBodyDynamics/left_foot/cartesianEndEffectorWrench:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.wL_WBDT[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S34>/left_foot_wrench

    /* Start for S-Function (WBToolbox): '<S82>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "imu_frame",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_k[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S82>/S-Function

    /* Start for S-Function (WBToolbox): '<S79>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_ao[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S79>/S-Function

    /* Start for S-Function (WBToolbox): '<S83>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "root_link",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_nj[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S83>/S-Function

    /* Start for S-Function (WBToolbox): '<S34>/inertial' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            12.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/inertial",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.inertial_n[0]), {1, 12});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S34>/inertial

    /* Start for S-Function (WBToolbox): '<S78>/Neck Position' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/head/state:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.NeckPosition_p[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S78>/Neck Position

    /* Start for S-Function (WBToolbox): '<S106>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "com", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.Switch6[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_jj[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S106>/S-Function

    /* Start for S-Function (WBToolbox): '<S92>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "imu_frame",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_oz[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S92>/S-Function

    /* Start for S-Function (WBToolbox): '<S81>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_i[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S81>/S-Function

    /* Start for S-Function (WBToolbox): '<S93>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "root_link",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_ef[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S93>/S-Function

    /* Start for S-Function (WBToolbox): '<S80>/Neck Position' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/head/state:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.NeckPosition_k[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S80>/Neck Position

    /* Start for S-Function (WBToolbox): '<S107>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "com", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.Switch6_g[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_dn[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S107>/S-Function

    /* Start for S-Function (WBToolbox): '<S76>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.w_H_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_e0[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S76>/S-Function

    /* Start for S-Function (WBToolbox): '<S77>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.w_H_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_da[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S77>/S-Function

    /* Start for S-Function (WBToolbox): '<S75>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Velocity",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_nr[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S75>/S-Function

    /* Start for S-Function (WBToolbox): '<S34>/Minimum Jerk Trajectory Generator' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 8.0, 1.0, 1.0, "ResetOnSettlingTimeChange"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ReadExternalSettlingTime"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "ReadInitialValue"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "ComputeFirstDerivative"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ComputeSecondDerivative"));
        params.storeParameter<double>(
            2.0, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 3.0, 1.0, 1.0, "SettlingTime"));
        params.storeParameter<double>(
            0.01, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 2.0, 1.0, 1.0, "SampleTime"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "MinimumJerkTrajectoryGenerator",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0,
            static_cast<void*>(&torqueBalancingYoga_B.TmpSignalConversionAtMinimumJer[0]),
            {1, 29});

        // Outputs
        blockInfo->setOutputSignal(
            0,
            static_cast<void*>(&torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator[0]),
            {1, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::MinimumJerkTrajectoryGenerator());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[0] =
            static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S34>/Minimum Jerk Trajectory Generator

    /* End of Start for SubSystem: '<S4>/State Machine Yoga' */

    /* Start for Enabled SubSystem: '<S4>/Internal Coordinator' */
    /* Start for S-Function (WBToolbox): '<S44>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "imu_frame",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_jo[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S44>/S-Function

    /* Start for S-Function (WBToolbox): '<S45>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "root_link",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_l[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S45>/S-Function

    /* Start for S-Function (WBToolbox): '<S42>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_k0[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S42>/S-Function

    /* Start for S-Function (WBToolbox): '<S43>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_bdr[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S43>/S-Function

    /* Start for S-Function (WBToolbox): '<S41>/Neck Position' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/head/state:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.NeckPosition_c[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S41>/Neck Position

    /* Start for S-Function (WBToolbox): '<S32>/IMU measurements' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            12.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/inertial",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.IMUmeasurements[0]), {1, 12});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S32>/IMU measurements

    /* Start for S-Function (WBToolbox): '<S57>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.Switch6_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_p[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S57>/S-Function

    /* Start for S-Function (WBToolbox): '<S58>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.Switch6_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_i0[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S58>/S-Function

    /* Start for Constant: '<S32>/Constant1' */
    torqueBalancingYoga_B.Constant1_c[0] =
        torqueBalancingYoga_P.Config.LEFT_RIGHT_FOOT_IN_CONTACT[0];
    torqueBalancingYoga_B.Constant1_c[1] =
        torqueBalancingYoga_P.Config.LEFT_RIGHT_FOOT_IN_CONTACT[1];

    /* Start for S-Function (WBToolbox): '<S56>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Velocity",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_m[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S56>/S-Function

    /* Start for S-Function (WBToolbox): '<S63>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "com", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.Switch6_b[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.jointAngles_p[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_ak[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S63>/S-Function

    /* End of Start for SubSystem: '<S4>/Internal Coordinator' */

    /* Start for S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator1' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 8.0, 1.0, 1.0, "ResetOnSettlingTimeChange"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ReadExternalSettlingTime"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "ReadInitialValue"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "ComputeFirstDerivative"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ComputeSecondDerivative"));
        params.storeParameter<double>(
            0.01,
            wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 3.0, 1.0, 1.0, "SettlingTime"));
        params.storeParameter<double>(
            0.01, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 2.0, 1.0, 1.0, "SampleTime"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "MinimumJerkTrajectoryGenerator",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[9]), {1, 23});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[110]), {1, 1});

        // Outputs
        blockInfo->setOutputSignal(
            0,
            static_cast<void*>(&torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator1[0]),
            {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::MinimumJerkTrajectoryGenerator());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[0] =
            static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator1

    /* Start for S-Function (WBToolbox): '<S18>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "MassMatrix",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_b[0]), {29, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::MassMatrix());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S18>/S-Function

    /* Start for S-Function (WBToolbox): '<S17>/S-Function' incorporates:
     *  Constant: '<S17>/Constant'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "InverseDynamics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[58]), {1, 6});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[64]), {1, 23});

        blockInfo->setInputSignal(
            4, static_cast<void*>(&torqueBalancingYoga_P.Constant_Value_d[0]), {1, 6});

        blockInfo->setInputSignal(5, static_cast<void*>(&torqueBalancingYoga_B.Gain[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_f[0]), {1, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::InverseDynamics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S17>/S-Function

    /* Start for S-Function (WBToolbox): '<S16>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "CentroidalMomentum",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[58]), {1, 6});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[64]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_c[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::CentroidalMomentum());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S16>/S-Function

    /* Start for S-Function (WBToolbox): '<S122>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "imu_frame",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_bd[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S122>/S-Function

    /* Start for S-Function (WBToolbox): '<S119>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_br[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S119>/S-Function

    /* Start for S-Function (WBToolbox): '<S123>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "root_link",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_d2[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S123>/S-Function

    /* Start for S-Function (WBToolbox): '<S109>/inertial' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.inertial_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            12.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/inertial",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.inertial[0]), {1, 12});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.inertial_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S109>/inertial

    /* Start for S-Function (WBToolbox): '<S118>/Neck Position' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/head/state:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.NeckPosition[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S118>/Neck Position

    /* Start for S-Function (WBToolbox): '<S132>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "imu_frame",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_o[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S132>/S-Function

    /* Start for S-Function (WBToolbox): '<S121>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_a[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S121>/S-Function

    /* Start for S-Function (WBToolbox): '<S133>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "root_link",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_P.Constant7_Value_o[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_c2[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S133>/S-Function

    /* Start for S-Function (WBToolbox): '<S120>/Neck Position' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<int>(
            6.0, wbt::ParameterMetadata(wbt::ParameterType::INT, 3.0, 1.0, 1.0, "SignalSize"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ErrorOnMissingPort"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "Autoconnect"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "WaitData"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ReadTimestamp"));
        params.storeParameter<double>(
            0.5, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 8.0, 1.0, 1.0, "Timeout"));
        params.storeParameter<std::string>(
            "/icub/head/state:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 2.0, 1.0, 1.0, "PortName"));
        params.storeParameter<std::string>(
            "YarpRead",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.NeckPosition_m[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpRead());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S120>/Neck Position

    /* Start for S-Function (WBToolbox): '<S114>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.Switch[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_e[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S114>/S-Function

    /* Start for S-Function (WBToolbox): '<S115>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.Switch[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_cf[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S115>/S-Function

    /* Start for S-Function (WBToolbox): '<S112>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "CentroidalMomentum",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.Switch[0]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.Switch5[0]), {1, 23});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.nu_b_equivalent[0]), {1, 6});

        blockInfo->setInputSignal(3, static_cast<void*>(&torqueBalancingYoga_B.Sum[0]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_on[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::CentroidalMomentum());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S112>/S-Function

    /* Start for S-Function (WBToolbox): '<S25>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_oe[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S25>/S-Function

    /* Start for S-Function (WBToolbox): '<S26>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_j[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S26>/S-Function

    /* Start for S-Function (WBToolbox): '<S23>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_c4[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S23>/S-Function

    /* Start for S-Function (WBToolbox): '<S24>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_b3[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S24>/S-Function

    /* Start for S-Function (WBToolbox): '<S20>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "l_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "DotJNu",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[58]), {1, 6});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[64]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_ej[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::DotJNu());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S20>/S-Function

    /* Start for S-Function (WBToolbox): '<S21>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "r_sole", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "DotJNu",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[58]), {1, 6});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[64]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_h[0]), {1, 6});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::DotJNu());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S21>/S-Function

    /* Start for S-Function (WBToolbox): '<S28>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "com", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "ForwardKinematics",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_dk[0]), {4, 4});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::ForwardKinematics());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S28>/S-Function

    /* Start for S-Function (WBToolbox): '<S22>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "com", wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "Frame"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "Jacobian",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[117]), {4, 4});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[87]), {1, 23});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_n[0]), {6, 29});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::Jacobian());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S22>/S-Function

    /* Start for S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator2' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 8.0, 1.0, 1.0, "ResetOnSettlingTimeChange"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "ReadExternalSettlingTime"));
        params.storeParameter<bool>(
            0.0,
            wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "ReadInitialValue"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "ComputeFirstDerivative"));
        params.storeParameter<bool>(
            1.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "ComputeSecondDerivative"));
        params.storeParameter<double>(
            0.01,
            wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 3.0, 1.0, 1.0, "SettlingTime"));
        params.storeParameter<double>(
            0.01, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 2.0, 1.0, 1.0, "SampleTime"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "MinimumJerkTrajectoryGenerator",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[0]), {1, 3});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.MultiportSwitch1[110]), {1, 1});

        // Outputs
        blockInfo->setOutputSignal(
            0,
            static_cast<void*>(&torqueBalancingYoga_B.MinimumJerkTrajectoryGenerator2[0]),
            {1, 3});

        blockInfo->setOutputSignal(
            1,
            static_cast<void*>(&torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_d[0]),
            {1, 3});

        blockInfo->setOutputSignal(
            2,
            static_cast<void*>(&torqueBalancingYoga_B.MinimumJerkTrajectoryGenerato_b[0]),
            {1, 3});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::MinimumJerkTrajectoryGenerator());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[0] =
            static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator2

    /* Start for Enabled SubSystem: '<S142>/One Foot' */

    /* Start for S-Function (WBToolbox): '<S144>/QP Two Feet' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "ComputeObjVal"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "UseUb"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "UseLb"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 2.0, 1.0, 1.0, "UseLbA"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "StopWhenFails"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 3.0, 1.0, 1.0, "UseUbA"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "QpOases",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.HessianMatrixQP1Foot[0]), {6, 6});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.gradientQP1Foot[0]), {1, 6});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.ConstraintsMatrixQP1Foot[0]), {19, 6});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.bVectorConstraintsQp1Foot[0]), {1, 19});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.QPTwoFeet_o1_d[0]), {1, 6});

        blockInfo->setOutputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.QPTwoFeet_o2_h), {1, 1});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::QpOases());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S144>/QP Two Feet

    /* End of Start for SubSystem: '<S142>/One Foot' */

    /* Start for Enabled SubSystem: '<S142>/Two Feet' */

    /* Start for S-Function (WBToolbox): '<S146>/QP Two Feet' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 6.0, 1.0, 1.0, "ComputeObjVal"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 5.0, 1.0, 1.0, "UseUb"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 4.0, 1.0, 1.0, "UseLb"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 2.0, 1.0, 1.0, "UseLbA"));
        params.storeParameter<bool>(
            0.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 7.0, 1.0, 1.0, "StopWhenFails"));
        params.storeParameter<bool>(
            1.0, wbt::ParameterMetadata(wbt::ParameterType::BOOL, 3.0, 1.0, 1.0, "UseUbA"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "QpOases",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.HessianMatrixQP2Feet[0]), {12, 12});

        blockInfo->setInputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.gradientQP2Feet[0]), {1, 12});

        blockInfo->setInputSignal(
            2, static_cast<void*>(&torqueBalancingYoga_B.ConstraintsMatrixQP2Feet[0]), {38, 12});

        blockInfo->setInputSignal(
            3, static_cast<void*>(&torqueBalancingYoga_B.bVectorConstraintsQp2Feet[0]), {1, 38});

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.QPTwoFeet_o1[0]), {1, 12});

        blockInfo->setOutputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.QPTwoFeet_o2), {1, 1});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::QpOases());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S146>/QP Two Feet

    /* End of Start for SubSystem: '<S142>/Two Feet' */

    /* Start for S-Function (WBToolbox): '<S5>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<double>(
            0.0,
            wbt::ParameterMetadata(
                wbt::ParameterType::DOUBLE, 5.0, 1.0, 1.0, "TrajectoryReference"));

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Torque",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "CtrlType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "SetReferences",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs
        blockInfo->setInputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.Saturation[0]), {1, 23});

        // Outputs

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::SetReferences());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S5>/S-Function

    /* Start for Enabled SubSystem: '<Root>/emergency stop: joint limits' */

    /* Start for S-Function (WBToolbox): '<S157>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "ControlBoardPosition",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "LimitType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetLimits",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_o1[0]), {1, 23});

        blockInfo->setOutputSignal(
            1, static_cast<void*>(&torqueBalancingYoga_B.SFunction_o2[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetLimits());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S157>/S-Function

    /* End of Start for SubSystem: '<Root>/emergency stop: joint limits' */

    /* Start for Enabled SubSystem: '<S2>/Visualizer' */

    /* Start for S-Function (WBToolbox): '<S12>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Torque",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_bm[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S12>/S-Function

    /* Start for S-Function (WBToolbox): '<S11>/S-Function' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;

        {
            std::vector<double> valueVector;
            valueVector.reserve(3.0);
            valueVector.push_back(0.0);
            valueVector.push_back(0.0);
            valueVector.push_back(-9.81);
            params.storeParameter<double>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_DOUBLE, 2.0, 1.0, 3.0, "GravityVector"));
        }

        params.storeParameter<std::string>(
            "Joints Position",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "MeasuredType"));
        params.storeParameter<std::string>(
            "WBT",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "LocalName"));
        params.storeParameter<std::string>(
            "icub",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "RobotName"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "GetMeasurement",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "torqueBalancingYoga/Configuration",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "ConfBlockName"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(23.0);
            valueVector.push_back("torso_pitch");
            valueVector.push_back("torso_roll");
            valueVector.push_back("torso_yaw");
            valueVector.push_back("l_shoulder_pitch");
            valueVector.push_back("l_shoulder_roll");
            valueVector.push_back("l_shoulder_yaw");
            valueVector.push_back("l_elbow");
            valueVector.push_back("r_shoulder_pitch");
            valueVector.push_back("r_shoulder_roll");
            valueVector.push_back("r_shoulder_yaw");
            valueVector.push_back("r_elbow");
            valueVector.push_back("l_hip_pitch");
            valueVector.push_back("l_hip_roll");
            valueVector.push_back("l_hip_yaw");
            valueVector.push_back("l_knee");
            valueVector.push_back("l_ankle_pitch");
            valueVector.push_back("l_ankle_roll");
            valueVector.push_back("r_hip_pitch");
            valueVector.push_back("r_hip_roll");
            valueVector.push_back("r_hip_yaw");
            valueVector.push_back("r_knee");
            valueVector.push_back("r_ankle_pitch");
            valueVector.push_back("r_ankle_roll");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 23.0, "ControlledJoints"));
        }

        params.storeParameter<std::string>(
            "model.urdf",
            wbt::ParameterMetadata(wbt::ParameterType::STRUCT_STRING, 2.0, 1.0, 1.0, "UrdfFile"));

        {
            std::vector<std::string> valueVector;
            valueVector.reserve(5.0);
            valueVector.push_back("torso");
            valueVector.push_back("left_arm");
            valueVector.push_back("right_arm");
            valueVector.push_back("left_leg");
            valueVector.push_back("right_leg");
            params.storeParameter<std::string>(
                valueVector,
                wbt::ParameterMetadata(
                    wbt::ParameterType::STRUCT_CELL_STRING, 2.0, 1.0, 5.0, "ControlBoardsNames"));
        }

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(
            0, static_cast<void*>(&torqueBalancingYoga_B.SFunction_ol[0]), {1, 23});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::GetMeasurement());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S11>/S-Function

    /* End of Start for SubSystem: '<S2>/Visualizer' */

    /* Start for Enabled SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' */

    /* Start for S-Function (WBToolbox): '<S159>/Simulator Synchronizer' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<double>(
            0.01, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 2.0, 1.0, 1.0, "Period"));
        params.storeParameter<std::string>(
            "/WBT_synchronizer/rpc:o",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 4.0, 1.0, 1.0, "RpcPort"));
        params.storeParameter<std::string>(
            "SimulatorSynchronizer",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));
        params.storeParameter<std::string>(
            "/clock/rpc",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 3.0, 1.0, 1.0, "GazeboClockPort"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::SimulatorSynchronizer());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[0] =
            static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S159>/Simulator Synchronizer

    /* End of Start for SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' */

    /* Start for Enabled SubSystem: '<S8>/REAL_TIME_SYNC' */

    /* Start for S-Function (WBToolbox): '<S160>/Real Time Synchronizer' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[1] =
            static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<double>(
            0.01, wbt::ParameterMetadata(wbt::ParameterType::DOUBLE, 2.0, 1.0, 1.0, "Period"));
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "RealTimeSynchronizer",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::RealTimeSynchronizer());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[0] =
            static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S160>/Real Time Synchronizer

    /* End of Start for SubSystem: '<S8>/REAL_TIME_SYNC' */

    /* Start for S-Function (WBToolbox): '<S8>/Yarp Clock' */
    {
        // Create and store the CoderBlockInformation object
        wbt::CoderBlockInformation* blockInfo = new wbt::CoderBlockInformation();
        torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[1] = static_cast<void*>(blockInfo);

        // Initialize the parameters
        // -------------------------
        wbt::Parameters params;
        params.storeParameter<std::string>(
            "WBToolboxLibrary",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 1.0, 1.0, 1.0, "libName"));
        params.storeParameter<std::string>(
            "YarpClock",
            wbt::ParameterMetadata(wbt::ParameterType::STRING, 0.0, 1.0, 1.0, "className"));

        // Store the parameters in the CoderBlockInformation object
        blockInfo->storeRTWParameters(params);

        // Initialize input / output Signals
        // ---------------------------------

        // Inputs

        // Outputs
        blockInfo->setOutputSignal(0, static_cast<void*>(&torqueBalancingYoga_B.YarpClock), {1, 1});

        // Initialize the class
        // --------------------

        // Allocate the block object
        wbt::Block* blockPtr = static_cast<wbt::Block*>(new wbt::YarpClock());

        // Run a dummy configureSizeAndPorts step. This is currently required for properly
        // handling optional input / outputs static variables.
        // TODO: find a better way to handle them.
        {
            auto tmpCoderBlockInfo =
                std::unique_ptr<wbt::CoderBlockInformation>(new wbt::CoderBlockInformation);
            tmpCoderBlockInfo->storeRTWParameters(params);
            blockPtr->configureSizeAndPorts(tmpCoderBlockInfo.get());
        }

        // Initialize the block
        bool ok = blockPtr->initialize(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Initialize]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Call the initializeInitialConditions() method
        ok = blockPtr->initializeInitialConditions(blockInfo);

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[InitializeInitialConditions]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }

        // Store the block in the PWork vector
        torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[0] = static_cast<void*>(blockPtr);
    }

    // End of S-Function Block: <S8>/Yarp Clock

    /* SystemInitialize for Enabled SubSystem: '<S4>/State Machine Yoga' */
    /* SystemInitialize for MATLAB Function: '<S85>/MATLAB Function' */
    torqueBalan_MATLABFunction_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_i);

    /* SystemInitialize for MATLAB Function: '<S86>/MATLAB Function' */
    torqueBal_MATLABFunction_m_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_oa);

    /* SystemInitialize for MATLAB Function: '<S68>/MATLAB Function' */
    torqueBal_MATLABFunction_l_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_f);

    /* SystemInitialize for MATLAB Function: '<S69>/MATLAB Function' */
    torqueBa_MATLABFunction_mb_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_d);

    /* SystemInitialize for MATLAB Function: '<S95>/MATLAB Function' */
    torqueBalan_MATLABFunction_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_ad);

    /* SystemInitialize for MATLAB Function: '<S96>/MATLAB Function' */
    torqueBal_MATLABFunction_m_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_j);

    /* SystemInitialize for MATLAB Function: '<S34>/stateMachineYogaFCN' */
    torqueBalancingYoga_DW.state_not_empty = false;
    torqueBalancingYoga_DW.tSwitch_not_empty = false;
    torqueBalancingYoga_DW.w_H_fixedLink_not_empty = false;
    torqueBalancingYoga_DW.secondYoga_not_empty = false;

    /* End of SystemInitialize for SubSystem: '<S4>/State Machine Yoga' */

    /* SystemInitialize for Enabled SubSystem: '<S4>/Internal Coordinator' */

    /* SystemInitialize for MATLAB Function: '<S47>/MATLAB Function' */
    torqueBalan_MATLABFunction_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_os);

    /* SystemInitialize for MATLAB Function: '<S48>/MATLAB Function' */
    torqueBal_MATLABFunction_m_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_l5);

    /* SystemInitialize for MATLAB Function: '<S39>/MATLAB Function' */
    torqueBal_MATLABFunction_l_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_a);

    /* SystemInitialize for MATLAB Function: '<S38>/MATLAB Function' */
    torqueBa_MATLABFunction_mb_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_k);

    /* End of SystemInitialize for SubSystem: '<S4>/Internal Coordinator' */

    /* SystemInitialize for MATLAB Function: '<S125>/MATLAB Function' */
    torqueBalan_MATLABFunction_Init(&torqueBalancingYoga_DW.sf_MATLABFunction);

    /* SystemInitialize for MATLAB Function: '<S126>/MATLAB Function' */
    torqueBal_MATLABFunction_m_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_h);

    /* SystemInitialize for MATLAB Function: '<S135>/MATLAB Function' */
    torqueBalan_MATLABFunction_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_o);

    /* SystemInitialize for MATLAB Function: '<S136>/MATLAB Function' */
    torqueBal_MATLABFunction_m_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_l);

    /* SystemInitialize for MATLAB Function: '<S162>/MATLAB Function' */
    torqueBa_MATLABFunction_mb_Init(&torqueBalancingYoga_DW.sf_MATLABFunction_p);

    /* SystemInitialize for MATLAB Function: '<S9>/Saturate the Torque Derivative' */
    torqueBalancingYoga_DW.uPrev_not_empty = false;
}

/* Model terminate function */
void torqueBalancingYogaModelClass::terminate()
{
    /* Terminate for S-Function (WBToolbox): '<S156>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S156>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S31>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_c.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S31>/S-Function

    /* Terminate for Enabled SubSystem: '<S4>/State Machine Yoga' */

    /* Terminate for S-Function (WBToolbox): '<S34>/right_foot_wrench' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.right_foot_wrench_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S34>/right_foot_wrench

    /* Terminate for S-Function (WBToolbox): '<S34>/left_foot_wrench' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.left_foot_wrench_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S34>/left_foot_wrench

    /* Terminate for S-Function (WBToolbox): '<S82>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_pz.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S82>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S79>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bg.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S79>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S83>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_m.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S83>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S34>/inertial' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.inertial_PWORK_a.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S34>/inertial

    /* Terminate for S-Function (WBToolbox): '<S78>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_n.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S78>/Neck Position

    /* Terminate for S-Function (WBToolbox): '<S106>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hk.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S106>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S92>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_kk.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S92>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S81>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_lc.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S81>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S93>/S-Function' incorporates:
     *  Constant: '<S66>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_m2.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S93>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S80>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_p.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S80>/Neck Position

    /* Terminate for S-Function (WBToolbox): '<S107>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d0.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S107>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S76>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j2.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S76>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S77>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_jg.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S77>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S75>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bs.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S75>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S34>/Minimum Jerk Trajectory Generator' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator_.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S34>/Minimum Jerk Trajectory Generator

    /* End of Terminate for SubSystem: '<S4>/State Machine Yoga' */

    /* Terminate for Enabled SubSystem: '<S4>/Internal Coordinator' */

    /* Terminate for S-Function (WBToolbox): '<S44>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ga.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S44>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S45>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_es.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S45>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S42>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_o.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S42>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S43>/S-Function' incorporates:
     *  Constant: '<S35>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_kj.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S43>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S41>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_c.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S41>/Neck Position

    /* Terminate for S-Function (WBToolbox): '<S32>/IMU measurements' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.IMUmeasurements_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S32>/IMU measurements

    /* Terminate for S-Function (WBToolbox): '<S57>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_l5.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S57>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S58>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_i.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S58>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S56>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j1.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S56>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S63>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_gx.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S63>/S-Function

    /* End of Terminate for SubSystem: '<S4>/Internal Coordinator' */

    /* Terminate for S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator1' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator1.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator1

    /* Terminate for S-Function (WBToolbox): '<S18>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_n.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S18>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S17>/S-Function' incorporates:
     *  Constant: '<S17>/Constant'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_h.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S17>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S16>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_e.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S16>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S122>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_b.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S122>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S119>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S119>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S123>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_nq.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S123>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S109>/inertial' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.inertial_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.inertial_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S109>/inertial

    /* Terminate for S-Function (WBToolbox): '<S118>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S118>/Neck Position

    /* Terminate for S-Function (WBToolbox): '<S132>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_k.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S132>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S121>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_er.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S121>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S133>/S-Function' incorporates:
     *  Constant: '<S113>/Constant7'
     */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hg.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S133>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S120>/Neck Position' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.NeckPosition_PWORK_j.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S120>/Neck Position

    /* Terminate for S-Function (WBToolbox): '<S114>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_dc.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S114>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S115>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_a.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S115>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S112>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_c4.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S112>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S25>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_j.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S25>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S26>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_p.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S26>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S23>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ne.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S23>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S24>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_bp.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S24>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S20>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d4.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S20>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S21>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_d1.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S21>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S28>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_hl.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S28>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S22>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_ay.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S22>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S33>/Minimum Jerk Trajectory Generator2' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.MinimumJerkTrajectoryGenerator2.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S33>/Minimum Jerk Trajectory Generator2

    /* Terminate for Enabled SubSystem: '<S142>/One Foot' */

    /* Terminate for S-Function (WBToolbox): '<S144>/QP Two Feet' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.QPTwoFeet_PWORK_e.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S144>/QP Two Feet

    /* End of Terminate for SubSystem: '<S142>/One Foot' */

    /* Terminate for Enabled SubSystem: '<S142>/Two Feet' */

    /* Terminate for S-Function (WBToolbox): '<S146>/QP Two Feet' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.QPTwoFeet_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S146>/QP Two Feet

    /* End of Terminate for SubSystem: '<S142>/Two Feet' */

    /* Terminate for S-Function (WBToolbox): '<S5>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_l.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S5>/S-Function

    /* Terminate for Enabled SubSystem: '<Root>/emergency stop: joint limits' */

    /* Terminate for S-Function (WBToolbox): '<S157>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_g.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S157>/S-Function

    /* End of Terminate for SubSystem: '<Root>/emergency stop: joint limits' */

    /* Terminate for Enabled SubSystem: '<S2>/Visualizer' */

    /* Terminate for S-Function (WBToolbox): '<S12>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_f.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S12>/S-Function

    /* Terminate for S-Function (WBToolbox): '<S11>/S-Function' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr =
            static_cast<wbt::Block*>(torqueBalancingYoga_DW.SFunction_PWORK_dk.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S11>/S-Function

    /* End of Terminate for SubSystem: '<S2>/Visualizer' */

    /* Terminate for Enabled SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' */

    /* Terminate for S-Function (WBToolbox): '<S159>/Simulator Synchronizer' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.SimulatorSynchronizer_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S159>/Simulator Synchronizer

    /* End of Terminate for SubSystem: '<S8>/GAZEBO_SYNCHRONIZER' */

    /* Terminate for Enabled SubSystem: '<S8>/REAL_TIME_SYNC' */

    /* Terminate for S-Function (WBToolbox): '<S160>/Real Time Synchronizer' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(
            torqueBalancingYoga_DW.RealTimeSynchronizer_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S160>/Real Time Synchronizer

    /* End of Terminate for SubSystem: '<S8>/REAL_TIME_SYNC' */

    /* Terminate for S-Function (WBToolbox): '<S8>/Yarp Clock' */
    {
        // Get the CoderBlockInformation from the PWork
        wbt::CoderBlockInformation* blockInfo = nullptr;
        blockInfo = static_cast<wbt::CoderBlockInformation*>(
            torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[1]);

        // Get the Block from the PWork
        wbt::Block* blockPtr = nullptr;
        blockPtr = static_cast<wbt::Block*>(torqueBalancingYoga_DW.YarpClock_PWORK.blockPWork[0]);

        // Terminate the class
        // -------------------
        bool ok;
        ok = blockPtr->terminate(blockInfo);
        delete blockInfo;
        delete blockPtr;

        // Report errors
        if (!ok) {
            std::string error = wbt::Log::getSingleton().getErrors();
            error = "[Terminate]" + error;

            // Trim the message if needed
            if (error.length() >= 1024) {
                error = error.substr(0, 1024 - 1);
            }

            // This shouldn't happen
            if (getRTM()->errorStatus) {
                delete getRTM()->errorStatus;
                getRTM()->errorStatus = nullptr;
            }

            getRTM()->errorStatus = new char[1024];
            sprintf(const_cast<char_T*>(getRTM()->errorStatus), "%s", error.c_str());
            return;
        }
    }

    // End of S-Function Block: <S8>/Yarp Clock
}

/* Constructor */
torqueBalancingYogaModelClass::torqueBalancingYogaModelClass()
{
    static const P_torqueBalancingYoga_T torqueBalancingYoga_P_temp = {
        /* Variable: Sm
         * Referenced by:
         *   '<S32>/joints.smoothingTime'
         *   '<S34>/stateMachineYogaFCN'
         */
        {1.0,
         2.0,
         2.0,

         {1.0, 1.0, 1.0, 0.9, 2.0, 2.0, 1.0, 1.0, 1.0, 0.9, 2.0, 2.0, 5.0},
         0.9,
         5.0,
         50.0,
         100.0,
         0.02,
         15.0,
         50.0,
         1.0,

         {0.0, 0.0, 0.0,   0.0,   0.0, 0.02,  0.0, 0.0,   0.0,    0.0,    0.0, 0.02,  0.0,
          0.0, 0.0, 0.005, 0.005, 0.0, -0.09, 0.0, 0.012, -0.015, -0.017, 0.0, 0.025, 0.0,
          0.0, 0.0, 0.0,   0.0,   0.0, 0.0,   0.0, 0.0,   0.0,    0.0,    0.0, 0.0,   0.0},
         1.0,
         1.0,
         1,
         0,
         0,
         0,
         0,
         0,
         0,
         0.6,
         0.6,

         {0.0,     -0.0348, 0.0864,  0.0,     -0.0348, 0.0864,  0.0,     0.0864,  0.0864,  0.0,
          -0.0348, 0.0864,  0.0,     0.0,     0.0779,  0.0258,  0.0,     0.0779,  0.0258,  0.0,
          0.0258,  0.0258,  0.0,     0.0779,  0.0258,  0.0,     0.0,     0.0429,  0.0152,  0.0,
          0.0429,  0.0152,  0.0,     0.0152,  0.0152,  0.0,     0.0429,  0.0152,  0.0,     0.0,
          -0.1493, 0.1253,  0.0,     -0.1493, 0.1253,  0.0,     0.1253,  0.1253,  0.0,     -0.1493,
          0.1253,  0.0,     0.0,     0.858,   0.8135,  0.0,     0.858,   0.8135,  0.0,     0.8135,
          0.8135,  0.0,     0.858,   0.8135,  0.0,     0.0,     0.2437,  0.3051,  0.0,     0.2437,
          0.3051,  0.0,     0.3051,  0.3051,  0.0,     0.2437,  0.3051,  0.0,     0.0,     0.871,
          0.7928,  0.0,     0.871,   0.7928,  0.0,     0.7928,  0.7928,  0.0,     0.871,   0.7928,
          0.0,     0.0,     -0.1493, 0.0563,  0.0,     -0.1493, 0.0563,  0.0,     0.0563,  0.0563,
          0.0,     -0.1493, 0.0563,  0.0,     0.0,     0.858,   0.6789,  0.0,     0.858,   0.6789,
          0.0,     0.6789,  0.6789,  0.0,     0.858,   0.6789,  0.0,     0.0,     0.2437,  0.334,
          0.0,     0.2437,  0.334,   0.0,     0.334,   0.334,   0.0,     0.2437,  0.334,   0.0,
          0.0,     0.871,   0.6214,  0.0,     0.871,   0.6214,  0.0,     0.6214,  0.6214,  0.0,
          0.871,   0.6214,  0.0,     0.0,     -0.0015, -0.0015, 0.0,     -0.0015, 0.0107,  0.0,
          0.0107,  0.0005,  0.0,     0.0005,  -0.0026, 0.0,     0.0,     -0.1109, -0.1109, 0.0,
          -0.1109, -0.0741, 0.0,     -0.0741, 0.0793,  0.0,     0.0793,  0.0225,  0.0,     0.0,
          -0.0001, -0.0001, 0.0,     -0.0001, -0.0001, 0.0,     -0.0001, -0.0014, 0.0,     -0.0014,
          0.0093,  0.0,     0.0,     0.0003,  0.0003,  0.0,     0.0003,  -0.012,  0.0,     -0.012,
          -0.0051, 0.0,     -0.0051, -0.002,  0.0,     0.0,     0.016,   0.016,   0.0,     0.016,
          0.0252,  0.0,     0.0252,  0.0073,  0.0,     0.0073,  0.0027,  0.0,     0.0,     0.163,
          0.163,   0.0,     0.163,   0.1369,  0.0,     0.1369,  -0.1151, 0.0,     -0.1151, -0.0277,
          0.0,     0.0,     0.0005,  0.0005,  0.0,     0.0005,  -0.0026, 0.0,     -0.0026, -0.0015,
          0.0,     -0.0015, 0.0107,  0.0,     0.0,     0.0793,  0.0793,  0.0,     0.0793,  0.0225,
          0.0,     0.0225,  -0.1109, 0.0,     -0.1109, -0.0741, 0.0,     0.0,     -0.0014, -0.0014,
          0.0,     -0.0014, 0.0093,  0.0,     0.0093,  -0.0001, 0.0,     -0.0001, -0.0001, 0.0,
          0.0,     -0.0051, -0.0051, 0.0,     -0.0051, -0.002,  0.0,     -0.002,  0.0003,  0.0,
          0.0003,  -0.012,  0.0,     0.0,     0.0073,  -0.106,  0.0,     0.0073,  0.0027,  0.0,
          0.0027,  0.016,   0.0,     0.016,   0.0252,  0.0,     0.0,     -0.1151, -0.1151, 0.0,
          -0.1151, -0.0277, 0.0,     -0.0277, 0.163,   0.0,     0.163,   0.1369,  0.0},

         {0.0,     0.9,
          1.8,     2.7,
          3.6,     4.5,
          5.4,     6.3,
          7.2,     8.1,
          9.0,     9.9,
          10.8,    11.700000000000001,
          12.6,    13.5,
          14.4,    15.3,
          16.2,    17.1,
          18.0,    18.900000000000002,
          19.8,    20.7,
          21.6,    22.5,
          -0.079,  -0.079,
          -0.0852, -0.0852,
          -0.079,  -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          0.2279,  0.1279,
          -0.3273, -0.4273,
          -0.2273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          -0.4273, -0.4273,
          0.4519,  0.4519,
          0.0821,  0.0821,
          0.4519,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          0.0821,  0.0821,
          -1.1621, -1.1621,
          0.1391,  0.1391,
          -1.1621, 0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.6663,  0.6663,
          1.4585,  1.4585,
          0.6663,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          0.4919,  0.4965,
          0.2464,  0.2464,
          0.4965,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.9947,  0.9947,
          0.3042,  0.3042,
          0.9947,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          -1.0717, -1.0717,
          -0.4181, -0.4181,
          -1.0717, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          1.2904,  1.2904,
          1.68,    1.68,
          1.2904,  1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          -0.2447, -0.2493,
          0.7373,  0.7373,
          -0.2493, 0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          1.0948,  1.0948,
          0.3031,  0.3031,
          1.0948,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.206,   0.206,
          0.206,   0.3473,
          0.4473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.3484,  0.3714,
          0.3714,  0.3514,
          0.3514,  0.3514,
          0.3514,  0.3514,
          0.3514,  0.3514,
          0.8514,  0.8514,
          0.8514,  0.8514,
          1.5514,  0.2514,
          -0.3514, 0.3514,
          0.8514,  0.8514,
          0.8514,  0.8514,
          1.5514,  0.2514,
          -0.3514, 0.3514,
          0.4008,  0.9599,
          0.9599,  1.3107,
          1.3107,  1.3107,
          1.3107,  1.3107,
          0.0107,  1.3107,
          1.3107,  0.3107,
          1.3107,  0.0107,
          0.3107,  0.0107,
          0.3107,  1.3107,
          1.3107,  0.3107,
          1.3107,  0.0107,
          0.3107,  0.0107,
          0.3107,  1.3107,
          -0.0004, 1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          -0.3672, -1.6594,
          -1.6594, -0.0189,
          -0.0189, -0.0189,
          -1.6217, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.106,  -0.106,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          -0.0875, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614},

         {0.0,     0.9,
          1.8,     2.7,
          3.6,     4.5,
          5.4,     6.3,
          7.2,     8.1,
          9.0,     9.9,
          10.8,    11.700000000000001,
          12.6,    13.5,
          14.4,    15.3,
          16.2,    17.1,
          18.0,    18.900000000000002,
          19.8,    20.7,
          21.6,    22.5,
          -0.079,  -0.079,
          -0.0852, -0.0852,
          -0.079,  -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.0852, -0.0852,
          -0.2279, -0.1279,
          0.3273,  0.4273,
          0.2273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          0.4273,  0.4273,
          -0.4519, -0.4519,
          -0.0821, -0.0821,
          -0.4519, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -0.0821, -0.0821,
          -1.0717, -1.0717,
          -0.4181, -0.4181,
          -1.0717, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          -0.4181, -0.4181,
          1.2904,  1.2904,
          1.68,    1.68,
          1.2904,  1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          1.68,    1.68,
          -0.2447, -0.2493,
          0.7373,  0.7373,
          -0.2493, 0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          0.7373,  0.7373,
          1.0948,  1.0948,
          0.3031,  0.3031,
          1.0948,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          0.3031,  0.3031,
          -1.1621, -1.1621,
          0.1391,  0.1391,
          -1.1621, 0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.1391,  0.1391,
          0.6663,  0.6663,
          1.4585,  1.4585,
          0.6663,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          1.4585,  1.4585,
          0.4919,  0.4965,
          0.2464,  0.2464,
          0.4965,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.2464,  0.2464,
          0.9947,  0.9947,
          0.3042,  0.3042,
          0.9947,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3042,  0.3042,
          0.3484,  0.3714,
          0.3714,  0.3514,
          0.3514,  0.3514,
          0.3514,  0.3514,
          0.3514,  0.3514,
          0.8514,  0.8514,
          0.8514,  0.8514,
          1.5514,  0.2514,
          -0.3514, 0.3514,
          0.8514,  0.8514,
          0.8514,  0.8514,
          1.5514,  0.2514,
          -0.3514, 0.3514,
          0.4008,  0.9599,
          0.9599,  1.3107,
          1.3107,  1.3107,
          1.3107,  1.3107,
          0.0107,  1.3107,
          1.3107,  0.3107,
          1.3107,  0.0107,
          0.3107,  0.0107,
          0.3107,  1.3107,
          1.3107,  0.3107,
          1.3107,  0.0107,
          0.3107,  0.0107,
          0.3107,  1.3107,
          -0.0004, 1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          1.3253,  1.3253,
          -0.3672, -1.6594,
          -1.6594, -0.0189,
          -0.0189, -0.0189,
          -1.6217, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.0189, -0.0189,
          -0.106,  -0.106,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          0.5,     0.5,
          -0.0875, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          -0.0614, -0.0614,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.2092,  0.2092,
          0.206,   0.206,
          0.206,   0.3473,
          0.4473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.6473,  0.6473,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          0.0006,  0.0006,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1741, -0.1741,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          -0.1044, -0.1044,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07,
          0.07,    0.07},

         {0.9,
          0.6,
          1.2,
          1.7999999999999998,
          2.4,
          3.0,
          3.5999999999999996,
          4.2,
          4.8,
          5.3999999999999995,
          6.0,
          6.6,
          7.1999999999999993,
          7.8,
          8.4,
          9.0,
          9.6,
          10.2,
          10.799999999999999,
          11.4,
          12.0,
          12.6,
          13.2,
          13.799999999999999,
          14.399999999999999,
          15.0},

         {0.9,
          0.6,
          1.2,
          1.7999999999999998,
          2.4,
          3.0,
          3.5999999999999996,
          4.2,
          4.8,
          5.3999999999999995,
          6.0,
          6.6,
          7.1999999999999993,
          7.8,
          8.4,
          9.0,
          9.6,
          10.2,
          10.799999999999999,
          11.4,
          12.0,
          12.6,
          13.2,
          13.799999999999999,
          14.399999999999999,
          15.0}},

        /* Variable: Config
         * Referenced by:
         *   '<Root>/ON_GAZEBO 1'
         *   '<S2>/ON_GAZEBO 3'
         *   '<S2>/ON_GAZEBO 4'
         *   '<S8>/ON_GAZEBO '
         *   '<S9>/Saturate the Torque Derivative'
         *   '<S9>/Constant'
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
         *   '<S111>/Constant'
         *   '<S111>/Gain'
         *   '<S15>/Add motor reflected inertias'
         *   '<S35>/Constant1'
         *   '<S142>/ON_GAZEBO 3'
         *   '<S142>/ON_GAZEBO 4'
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
         *   '<S144>/Process QP output'
         *   '<S146>/Process QP output'
         *   '<S49>/USE_IMU4EST_BASE1'
         *   '<S87>/USE_IMU4EST_BASE1'
         *   '<S97>/USE_IMU4EST_BASE1'
         *   '<S127>/USE_IMU4EST_BASE1'
         *   '<S137>/USE_IMU4EST_BASE1'
         */
        {600.0,
         1,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0.01,
         0,

         {3.0, 4.0, 4.0, 6.0, 6.0},
         1,
         1,
         1,
         0,
         1,
         1,
         1,
         1,

         {1.0, 1.0},
         1,
         0,
         1,
         1,
         300.0,
         0.0,

         {0.0, 0.0, 0.0},
         0.0,
         0.0,

         {0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.01},

         {8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          5.8500000000000007e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          8.27e-6},

         {0.0,
          0.0,
          0.54999999999999993,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          -0.5,
          0.5,
          0.27499999999999997,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.5,
          0.5,
          0.27499999999999997,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          -1.0,
          -1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          -0.615,
          0.615,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          -0.615,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.615,
          -0.615,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.615,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          1.0},
         0.0},

        /* Variable: Gain
         * Referenced by:
         *   '<S6>/Balancing Controller '
         *   '<S34>/stateMachineYogaFCN'
         */
        {{50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,  50.0,
          100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 150.0, 100.0, 100.0, 100.0, 100.0, 100.0,
          5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0},

         {0.94280904158206336, 0.94280904158206336, 0.94280904158206336, 0.94280904158206336,
          0.94280904158206336, 0.94280904158206336, 0.94280904158206336, 0.94280904158206336,
          0.94280904158206336, 0.94280904158206336, 0.94280904158206336, 0.94280904158206336,
          0.94280904158206336, 1.3333333333333333,  1.3333333333333333,  1.3333333333333333,
          1.3333333333333333,  1.3333333333333333,  1.3333333333333333,  1.3333333333333333,
          1.6329931618554521,  1.3333333333333333,  1.3333333333333333,  1.3333333333333333,
          1.3333333333333333,  1.3333333333333333,  0.29814239699997197, 0.29814239699997197,
          0.29814239699997197, 0.29814239699997197, 0.29814239699997197, 0.29814239699997197,
          0.29814239699997197, 0.29814239699997197, 0.29814239699997197, 0.29814239699997197,
          0.29814239699997197, 0.29814239699997197, 0.29814239699997197},
         3.0,
         0.69282032302755092,

         {10.0,  10.0,  10.0,  90.0,  30.0,  30.0,  10.0,  10.0,  10.0,  90.0,  30.0,  30.0,
          90.0,  30.0,  30.0,  30.0,  90.0,  30.0,  30.0,  30.0,  30.0,  30.0,  90.0,  30.0,
          30.0,  120.0, 20.0,  20.0,  20.0,  90.0,  30.0,  30.0,  20.0,  20.0,  20.0,  90.0,
          30.0,  30.0,  90.0,  10.0,  10.0,  10.0,  10.0,  5.0,   10.0,  10.0,  10.0,  10.0,
          10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  5.0,   10.0,  10.0,  10.0,
          10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  20.0,  10.0,
          10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  8.0,   8.0,   8.0,   10.0,  10.0,  10.0,
          8.0,   8.0,   8.0,   10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,
          10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,
          10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,
          10.0,  20.0,  20.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  8.0,   8.0,
          8.0,   10.0,  10.0,  10.0,  8.0,   8.0,   8.0,   10.0,  10.0,  10.0,  10.0,  30.0,
          30.0,  30.0,  100.0, 200.0, 100.0, 30.0,  30.0,  30.0,  100.0, 220.0, 220.0, 220.0,
          30.0,  30.0,  50.0,  200.0, 250.0, 350.0, 50.0,  50.0,  50.0,  50.0,  550.0, 550.0,
          550.0, 20.0,  20.0,  30.0,  100.0, 20.0,  20.0,  60.0,  60.0,  30.0,  30.0,  220.0,
          220.0, 220.0, 20.0,  20.0,  60.0,  400.0, 20.0,  200.0, 30.0,  30.0,  60.0,  50.0,
          200.0, 200.0, 200.0, 100.0, 100.0, 100.0, 100.0, 10.0,  10.0,  5.0,   100.0, 100.0,
          100.0, 65.0,  65.0,  65.0,  100.0, 100.0, 100.0, 100.0, 10.0,  100.0, 5.0,   100.0,
          100.0, 100.0, 300.0, 300.0, 300.0, 30.0,  30.0,  30.0,  100.0, 220.0, 440.0, 30.0,
          30.0,  30.0,  100.0, 200.0, 100.0, 100.0, 50.0,  50.0,  30.0,  50.0,  550.0, 1100.0,
          30.0,  30.0,  30.0,  200.0, 250.0, 350.0, 350.0, 30.0,  30.0,  20.0,  30.0,  220.0,
          440.0, 30.0,  30.0,  20.0,  100.0, 20.0,  20.0,  20.0,  60.0,  60.0,  20.0,  50.0,
          200.0, 400.0, 20.0,  20.0,  20.0,  100.0, 20.0,  200.0, 200.0, 100.0, 100.0, 100.0,
          100.0, 65.0,  130.0, 5.0,   100.0, 100.0, 10.0,  10.0,  10.0,  10.0,  100.0, 100.0,
          100.0, 100.0, 300.0, 600.0, 5.0,   100.0, 100.0, 10.0,  10.0,  100.0, 100.0},

         {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
         2.0},

        /* Variable: Reg
         * Referenced by:
         *   '<S6>/Balancing Controller '
         *   '<S109>/References for L'
         *   '<S36>/Compute Base Velocity'
         *   '<S65>/Compute Base Velocity'
         */
        {1.0e-7, 0.07, 1.0e-5, 0.1, 0.0, 1.0e-7},

        /* Variable: ConstraintsMatrix
         * Referenced by: '<S6>/Constant'
         */
        {3.7320508075688785,
         1.0000000000000002,
         0.26794919243112281,
         -0.2679491924311227,
         -0.99999999999999911,
         -3.7320508075688745,
         -3.7320508075688807,
         -1.0000000000000011,
         -0.26794919243112292,
         0.2679491924311222,
         0.99999999999999878,
         3.7320508075688732,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         1.0,
         1.0,
         1.0,
         1.0,
         1.0,
         1.0,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         -1.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -1.2440169358562927,
         -0.45534180126147961,
         -0.33333333333333331,
         -0.33333333333333331,
         -0.45534180126147927,
         -1.2440169358562914,
         -1.2440169358562936,
         -0.45534180126147983,
         -0.33333333333333331,
         -0.33333333333333331,
         -0.45534180126147922,
         -1.2440169358562909,
         -0.013333333333333334,
         -0.013333333333333334,
         -1.0,
         -0.07,
         -0.12,
         -0.045,
         -0.05,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -1.0,
         1.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         0.0,
         0.0,
         0.0,
         1.0,
         -1.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         -0.0,
         1.0,
         -1.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0},

        /* Variable: ROBOT_DOF_FOR_SIMULINK
         * Referenced by: '<S6>/    5'
         */
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},

        /* Variable: bVectorConstraints
         * Referenced by: '<S6>/Constant1'
         */
        {0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         -10.0,
         0.0,
         0.0,
         0.0,
         0.0},

        /* Variable: Sat
         * Referenced by: '<Root>/Saturation'
         */
        {60.0},

        /* Mask Parameter: CompareToConstant_const
         * Referenced by: '<S50>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_j
         * Referenced by: '<S52>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_o
         * Referenced by: '<S61>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_ok
         * Referenced by: '<S59>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_p
         * Referenced by: '<S88>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_od
         * Referenced by: '<S90>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_ph
         * Referenced by: '<S102>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_pf
         * Referenced by: '<S104>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_f
         * Referenced by: '<S98>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_e
         * Referenced by: '<S100>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_h
         * Referenced by: '<S128>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_c
         * Referenced by: '<S130>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_ea
         * Referenced by: '<S138>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_e5
         * Referenced by: '<S140>/Constant'
         */
        0.0,

        /* Mask Parameter: CompareToConstant_const_l
         * Referenced by: '<S163>/Constant'
         */
        0.0,

        /* Mask Parameter: Coordinator_BitMask
         * Referenced by: '<S4>/Coordinator'
         */
        1U,

        /* Mask Parameter: Yoga_BitMask
         * Referenced by: '<S4>/Yoga'
         */
        2U,

        /* Expression: pi/180
         * Referenced by: '<S49>/Gain'
         */
        0.017453292519943295,

        /* Expression: eye(4)
         * Referenced by: '<S35>/Constant7'
         */
        {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},

        /* Expression: zeros(3,1)
         * Referenced by: '<S49>/Constant'
         */
        {0.0, 0.0, 0.0},

        /* Expression: zeros(6,1)
         * Referenced by: '<S32>/Constant2'
         */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},

        /* Expression: Gain.impedances(1,:)
         * Referenced by: '<S32>/Constant3'
         */
        {10.0, 30.0, 20.0, 10.0,  10.0,  10.0, 8.0,  10.0, 10.0, 10.0,  8.0,  30.0,
         30.0, 20.0, 20.0, 100.0, 100.0, 30.0, 50.0, 30.0, 60.0, 100.0, 100.0},

        /* Expression: 1
         * Referenced by: '<S32>/Constant4'
         */
        1.0,

        /* Expression: diag(Gain.KP_COM)
         * Referenced by: '<S32>/Constant5'
         */
        {50.0, 100.0, 5.0},

        /* Expression: diag(Gain.KD_COM)
         * Referenced by: '<S32>/Constant6'
         */
        {0.94280904158206336, 1.3333333333333333, 0.29814239699997197},

        /* Expression: pi/180
         * Referenced by: '<S87>/Gain'
         */
        0.017453292519943295,

        /* Expression: pi/180
         * Referenced by: '<S97>/Gain'
         */
        0.017453292519943295,

        /* Expression: eye(4)
         * Referenced by: '<S66>/Constant7'
         */
        {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},

        /* Expression: zeros(3,1)
         * Referenced by: '<S87>/Constant'
         */
        {0.0, 0.0, 0.0},

        /* Expression: zeros(3,1)
         * Referenced by: '<S97>/Constant'
         */
        {0.0, 0.0, 0.0},

        /* Expression: zeros(6,1)
         * Referenced by: '<S34>/Constant1'
         */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},

        /* Expression: pi/180
         * Referenced by: '<S127>/Gain'
         */
        0.017453292519943295,

        /* Expression: pi/180
         * Referenced by: '<S137>/Gain'
         */
        0.017453292519943295,

        /* Expression: 0.01
         * Referenced by: '<S7>/index1'
         */
        0.01,

        /* Expression: zeros(6,1)
         * Referenced by: '<S17>/Constant'
         */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},

        /* Expression: 0
         * Referenced by: '<S17>/Gain'
         */
        0.0,

        /* Expression: eye(4)
         * Referenced by: '<S113>/Constant7'
         */
        {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},

        /* Expression: zeros(3,1)
         * Referenced by: '<S127>/Constant'
         */
        {0.0, 0.0, 0.0},

        /* Expression: zeros(3,1)
         * Referenced by: '<S137>/Constant'
         */
        {0.0, 0.0, 0.0},

        /* Expression: 0.1
         * Referenced by: '<S109>/Switch'
         */
        0.1,

        /* Computed Parameter: Constant_Value_n
         * Referenced by: '<S30>/Constant'
         */
        0U,

        /* Expression: Sm.SM_TYPE_BIN
         * Referenced by: '<S4>/Constant'
         */
        2U,

        /* Computed Parameter: Constant_Value_dx
         * Referenced by: '<S29>/Constant'
         */
        0U,

        /* Expression: Sm.SM_TYPE_BIN
         * Referenced by: '<S4>/Constant2'
         */
        2U}; /* Modifiable parameters */

    /* Initialize tunable parameters */
    torqueBalancingYoga_P = torqueBalancingYoga_P_temp;
}

/* Destructor */
torqueBalancingYogaModelClass::~torqueBalancingYogaModelClass()
{
    /* Currently there is no destructor body generated.*/
}

/* Real-Time Model get method */
RT_MODEL_torqueBalancingYoga_T* torqueBalancingYogaModelClass::getRTM()
{
    return (&torqueBalancingYoga_M);
}
