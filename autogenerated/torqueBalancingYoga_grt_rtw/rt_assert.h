/*
 * rt_assert.h
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
