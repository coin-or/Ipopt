// Copyright (C) 2005, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter    IBM       2005-08-12

#include "IpoptConfig.h"
#include "IpUtils.hpp"

#include <cstdlib>
#include <cmath>
#include <cfloat>
#ifdef HAVE_IEEFP_H
#include <ieeefp.h>
#endif
#include <ctime>
#include <cstdio>
#include <cstdarg>
#include <limits>

// The special treatment of vsnprintf on SUN has been suggested by Lou Hafer 2010/07/04
#if defined(HAVE_VSNPRINTF) && defined(__SUNPRO_CC)
namespace std
{
#include <iso/stdio_c99.h>
}
#endif

// The following code has been copied from CoinUtils' CoinTime

/** 8< (BEGIN) ******************************** */

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#else
// MacOS-X and FreeBSD needs sys/time.h
#if defined(__MACH__) || defined (__FreeBSD__)
#include <sys/time.h>
#endif
#if !defined(__MSVCRT__)
#include <sys/resource.h>
#endif
#endif

//#############################################################################

#if defined(_MSC_VER)

#if 0 // change this to 1 if want to use the win32 API
#include <windows.h>
#ifdef small
/* for some unfathomable reason (to me) rpcndr.h (pulled in by windows.h) does a
   '#define small char' */
#undef small
#endif
#define TWO_TO_THE_THIRTYTWO 4294967296.0
#define DELTA_EPOCH_IN_SECS  11644473600.0
inline double IpCoinGetTimeOfDay()
{
   FILETIME ft;

   GetSystemTimeAsFileTime(&ft);
   double t = ft.dwHighDateTime * TWO_TO_THE_THIRTYTWO + ft.dwLowDateTime;
   t = t / 10000000.0 - DELTA_EPOCH_IN_SECS;
   return t;
}
#else
#include <sys/types.h>
#include <sys/timeb.h>
inline double IpCoinGetTimeOfDay()
{
   struct _timeb timebuffer;
#pragma warning(disable:4996)
   _ftime( &timebuffer ); // C4996
#pragma warning(default:4996)
   return timebuffer.time + timebuffer.millitm / 1000.0;
}
#endif

#else

#include <sys/time.h>

inline double IpCoinGetTimeOfDay()
{
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec) / 1000000.0;
}

#endif // _MSC_VER

/** 8< (END) ******************************** */

namespace Ipopt
{

bool IsFiniteNumber(
   Number val
)
{
#ifdef IPOPT_C_FINITE
   return (bool)IPOPT_C_FINITE(val);
#else
   return true;
#endif

}

Number IpRandom01()
{
#ifdef IPOPT_HAS_DRAND48
   return Number(drand48());
#else
# ifdef IPOPT_HAS_RAND
   return Number(rand()) / Number(RAND_MAX);
# else
#  ifdef HAVE_STD__RAND
   return Number(std::rand()) / Number(RAND_MAX);
#  else
#   error "don't have function for random number generator"
#  endif
# endif
#endif
}

void IpResetRandom01()
{
#ifdef IPOPT_HAS_DRAND48
   srand48(1);
#else
# ifdef IPOPT_HAS_RAND
   srand(1);
# else
#  ifdef HAVE_STD__RAND
   std::srand(1);
#  else
#   error "don't have function for random number generator"
#  endif
# endif
#endif
}

static double Wallclock_firstCall_ = -1.;

// The following function were taken from CoinTime.hpp in COIN/Coin
Number CpuTime()
{
   double cpu_temp;

#if defined(_MSC_VER) || defined(__MSVCRT__)
   unsigned int ticksnow;        /* clock_t is same as int */

   ticksnow = (unsigned int)clock();

   cpu_temp = (double)((double)ticksnow / CLOCKS_PER_SEC);
#else
   struct rusage usage;
   getrusage(RUSAGE_SELF, &usage);
   cpu_temp = (double)usage.ru_utime.tv_sec;
   cpu_temp += 1.0e-6 * ((double) usage.ru_utime.tv_usec);
#endif

   return cpu_temp;
}

Number SysTime()
{
   double sys_temp;

#if defined(_MSC_VER) || defined(__MSVCRT__)
   // not yet implemented for Windows
   sys_temp = 0.;
#else
   struct rusage usage;
   getrusage(RUSAGE_SELF, &usage);
   sys_temp = (double)usage.ru_stime.tv_sec;
   sys_temp += 1.0e-6 * ((double) usage.ru_stime.tv_usec);
#endif

   return sys_temp;
}

Number WallclockTime()
{
   double callTime = IpCoinGetTimeOfDay();
   if( Wallclock_firstCall_ == -1. )
   {
      Wallclock_firstCall_ = callTime;
   }
   return callTime - Wallclock_firstCall_;
}

bool Compare_le(
   Number lhs,
   Number rhs,
   Number BasVal
)
{
   Number mach_eps = std::numeric_limits<Number>::epsilon();
   return (lhs - rhs <= 10.*mach_eps * std::abs(BasVal));
}

int Snprintf(
   char*       str,
   long        size,
   const char* format,
   ...
)
{
#if defined(HAVE_VSNPRINTF) && defined(__SUNPRO_CC)
   std::va_list ap;
#else
   va_list ap;
#endif
   va_start(ap, format);
   int ret;
#ifdef IPOPT_HAS_VA_COPY
   va_list apcopy;
   va_copy(apcopy, ap);
# ifdef HAVE_VSNPRINTF
#  ifdef __SUNPRO_CC
   ret = std::vsnprintf(str, size, format, apcopy);
#  else
   ret = vsnprintf(str, size, format, apcopy);
#  endif
# else
#  ifdef HAVE__VSNPRINTF
   ret = _vsnprintf(str, size, format, apcopy);
#  else
   ret = vsprintf(str, format, apcopy);
   (void) size;
#  endif
# endif
   va_end(apcopy);
#else
# ifdef HAVE_VSNPRINTF
#  ifdef __SUNPRO_CC
   ret = std::vsnprintf(str, size, format, ap);
#  else
   ret = vsnprintf(str, size, format, ap);
#  endif
# else
#  ifdef HAVE__VSNPRINTF
   ret = _vsnprintf(str, size, format, ap);
#  else
   ret = vsprintf(str, format, ap);
   (void) size;
#  endif
# endif
#endif
   va_end(ap);
   return ret;
}

} //namespace Ipopt
