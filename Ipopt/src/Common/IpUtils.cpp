// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter    IBM       2005-08-12

#include "IpoptConfig.h"
#include "IpUtils.hpp"

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif

#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif

#ifdef HAVE_CTIME
# include <ctime>
#else
# ifdef HAVE_TIME_H
#  include <time.h>
# else
#  error "don't have header file for time"
# endif
#endif

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#ifdef HAVE_CSTDARG
# include <cstdarg>
#else
# ifdef HAVE_STDARG_H
#  include <stdarg.h>
# else
#  error "don't have header file for stdarg"
# endif
#endif

#include <limits>

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
  t = t/10000000.0 - DELTA_EPOCH_IN_SECS;
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
  return timebuffer.time + timebuffer.millitm/1000.0;
}
#endif

#else

#include <sys/time.h>

inline double IpCoinGetTimeOfDay()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec/1000000.0;
}

#endif // _MSC_VER

/** 8< (END) ******************************** */

namespace Ipopt
{

  bool IsFiniteNumber(Number val)
  {
#ifdef MY_C_FINITE
    return (bool)MY_C_FINITE(val);
#else

    return true;
#endif

  }

  Number IpRandom01()
  {
#ifdef HAVE_DRAND48
    return Number(drand48());
#else
# ifdef HAVE_RAND
    return Number(rand())/Number(RAND_MAX);
# else
#  ifdef HAVE_STD__RAND
    return Number(std::rand())/Number(RAND_MAX);
#  else
#   error "don't have function for random number generator"
#  endif
# endif
#endif
  }

  void IpResetRandom01()
  {
#ifdef HAVE_DRAND48
    srand48(1);
#else
# ifdef HAVE_RAND
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

    cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else

    struct rusage usage;
    getrusage(RUSAGE_SELF,&usage);
    cpu_temp = (double)usage.ru_utime.tv_sec;
    cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
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
    getrusage(RUSAGE_SELF,&usage);
    sys_temp = (double)usage.ru_stime.tv_sec;
    sys_temp += 1.0e-6*((double) usage.ru_stime.tv_usec);
#endif

    return sys_temp;
  }

  double WallclockTime()
  {
    double callTime = IpCoinGetTimeOfDay();
    if (Wallclock_firstCall_ == -1.) {
      Wallclock_firstCall_ = callTime;
    }
    return callTime - Wallclock_firstCall_;
  }

  bool Compare_le(Number lhs, Number rhs, Number BasVal)
  {
    Number mach_eps = std::numeric_limits<Number>::epsilon();
    return (lhs - rhs <= 10.*mach_eps*fabs(BasVal));
  }

  int Snprintf(char* str, long size, const char* format, ...)
  {
    va_list ap;
    va_start(ap, format);
    int ret;
#ifdef HAVE_VA_COPY
    va_list apcopy;
    va_copy(apcopy, ap);
# ifdef HAVE_VSNPRINTF
    ret = vsnprintf(str, size, format, apcopy);
# else
#  ifdef HAVE__VSNPRINTF
    ret = _vsnprintf(str, size, format, apcopy);
#  else
    ret = vsprintf(str, format, apcopy);
#  endif
    va_end(apcopy);
# endif
#else
# ifdef HAVE_VSNPRINTF
    ret = vsnprintf(str, size, format, ap);
# else
#  ifdef HAVE__VSNPRINTF
    ret = _vsnprintf(str, size, format, ap);
#  else
    ret = vsprintf(str, format, ap);
#  endif
# endif
#endif
    va_end(ap);
    return ret;
  }

} //namespace Ipopt
