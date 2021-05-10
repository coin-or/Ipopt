// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPUTILS_HPP__
#define __IPUTILS_HPP__

// Standard Ip Include Files
#include "IpTypes.hpp"
#include "IpDebug.hpp"

namespace Ipopt
{

template<typename T>
inline T Max(
   T a,
   T b
)
{
   return ((a) > (b) ? (a) : (b));
}

template<typename T>
inline T Max(
   T a,
   T b,
   T c
)
{
   T max = Max(a, b);
   max = Max(max, c);
   return max;
}

template<typename T>
inline T Max(
   T a,
   T b,
   T c,
   T d
)
{
   T max = Max(a, b, c);
   return Max(max, d);
}

template<typename T>
inline T Min(
   T a,
   T b
)
{
   return ((a) < (b) ? (a) : (b));
}

template<typename T>
inline T Min(
   T a,
   T b,
   T c
)
{
   T min = Min(a, b);
   return Min(min, c);
}

template<typename T>
inline T Min(
   T a,
   T b,
   T c,
   T d
)
{
   T min = Min(a, b, c);
   return Min(min, d);
}

/** Function returning true iff the argument is a valid double number
 *  (not NaN or Inf). */
IPOPTLIB_EXPORT bool IsFiniteNumber(
   Number val
);

/** Function returning a random number between 0 and 1 */
IPOPTLIB_EXPORT Number IpRandom01();

/** Function resetting the random number generator */
IPOPTLIB_EXPORT void IpResetRandom01();

/** method determining CPU time */
IPOPTLIB_EXPORT Number CpuTime();

/** method determining system time */
IPOPTLIB_EXPORT Number SysTime();

/** method determining wallclock time since first call */
IPOPTLIB_EXPORT Number WallclockTime();

/** Method for comparing two numbers within machine precision.
 *
 *  @return true, if lhs is less or equal the rhs, relaxing
 *  this inequality by something a little larger than machine
 *  precision relative to the absolute value of BasVal
 */
IPOPTLIB_EXPORT bool Compare_le(
   Number lhs,
   Number rhs,
   Number BasVal
);

/** Method for printing a formatted output to a string with given size. */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
IPOPTLIB_EXPORT int Snprintf(
   char*       str,
   long        size,
   const char* format,
   ...
);

} //namespace Ipopt

#endif
