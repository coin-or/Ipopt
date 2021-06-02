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

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace Ipopt
{

template<typename T>
inline T Max(
   T a,
   T b
)
{
   return std::max(a, b);
}

template<typename T>
inline T Max(
   T a,
   T b,
   T c
)
{
   return std::max(std::max(a, b), c);
}

template<typename T>
inline T Max(
   T a,
   T b,
   T c,
   T d
)
{
   return std::max(std::max(a, b), std::max(c, d));
}

template<typename T>
inline T Min(
   T a,
   T b
)
{
   return std::min(a, b);
}

template<typename T>
inline T Min(
   T a,
   T b,
   T c
)
{
   return std::min(std::min(a, b), c);
}

template<typename T>
inline T Min(
   T a,
   T b,
   T c,
   T d
)
{
   return std::min(std::min(a, b), std::min(c, d));
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

/** Method to calculate new length for a memory increase based on a recommendation and limits in integer type.
 *
 * Checks whether recommended can be represented by T.
 * If so, sets len to min of recommendation and min.
 * If not, sets len to maximal value available for T, if this is larger than current value for len.
 * If not, throws a std::overflow_error exception.
 *
 * @since 3.14.0
 */
template<typename T>
inline void ComputeMemIncrease(
   T&          len,          ///< current length on input, new length on output
   double      recommended,  ///< recommended size
   T           min,          ///< minimal size that should ensured
   const char* context       ///< context from where this function is called - used to setup message for exception
)
{
   if( recommended >= std::numeric_limits<T>::max() )
   {
      // increase len to the maximum possible, if that is still an increase
      if( len < std::numeric_limits<T>::max() )
      {
         len = std::numeric_limits<T>::max();
      }
      else
      {
         DBG_ASSERT(context != NULL);
         std::stringstream what;
         what << "Cannot allocate more than " << std::numeric_limits<T>::max()*sizeof(T) << " bytes for " << context << " due to limitation on integer type";
         throw std::overflow_error(what.str());
      }
   }
   else
   {
      len = Max(min, (T) recommended);
   }
}

} //namespace Ipopt

#endif
