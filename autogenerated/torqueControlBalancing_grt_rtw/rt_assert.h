/*
 * rt_assert.h
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

#ifndef RTW_HEADER_rt_assert_h_
#define RTW_HEADER_rt_assert_h_

/*=========*
 * Asserts *
 *=========*/
#ifndef utAssert
#if defined(DOASSERTS)
#if !defined(PRINT_ASSERTS)
#include <assert.h>
#define utAssert(exp) assert(exp)
#else
#include <stdio.h>

static void _assert(char* statement, char* file, int line)
{
    printf("%s in %s on line %d\n", statement, file, line);
}

#define utAssert(_EX) ((_EX) ? (void) 0 : _assert(#_EX, __FILE__, __LINE__))
#endif

#else
#define utAssert(exp) /* do nothing */
#endif
#endif
#endif /* RTW_HEADER_rt_assert_h_ */
