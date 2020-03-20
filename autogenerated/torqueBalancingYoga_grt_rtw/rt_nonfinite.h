/*
 * rt_nonfinite.h
 *
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, educational organizations only. Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "torqueBalancingYoga".
 *
 * Model version              : 1.3291
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C++ source code generated on : Thu Sep 27 14:12:37 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_rt_nonfinite_h_
#define RTW_HEADER_rt_nonfinite_h_
#include "rtwtypes.h"
#include <stddef.h>
#ifdef __cplusplus

extern "C" {

#endif

extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
extern void rt_InitInfAndNaN(size_t realSize);
extern boolean_T rtIsInf(real_T value);
extern boolean_T rtIsInfF(real32_T value);
extern boolean_T rtIsNaN(real_T value);
extern boolean_T rtIsNaNF(real32_T value);
typedef struct
{
    struct
    {
        uint32_T wordH;
        uint32_T wordL;
    } words;
} BigEndianIEEEDouble;

typedef struct
{
    struct
    {
        uint32_T wordL;
        uint32_T wordH;
    } words;
} LittleEndianIEEEDouble;

typedef struct
{
    union {
        real32_T wordLreal;
        uint32_T wordLuint;
    } wordL;
} IEEESingle;

#ifdef __cplusplus

} /* extern "C" */
#endif
#endif /* RTW_HEADER_rt_nonfinite_h_ */
