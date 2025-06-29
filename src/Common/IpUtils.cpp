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
#include <csignal>
#include <limits>

#if defined(_MSC_VER) && _MSC_VER < 1900
#define vsnprintf _vsnprintf
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
#ifdef max
#undef max
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
#  ifdef IPOPT_HAS_STD__RAND
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
#  ifdef IPOPT_HAS_STD__RAND
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

static bool registered_handler = false;
static unsigned int abortcountdown_ = std::numeric_limits<unsigned int>::max();
static void (*handle_interrupt_)(void) = NULL;
static bool* interrupt_flag_ = NULL;

static void sighandler(
   int /* signum */
)
{
   if( interrupt_flag_ != NULL )
   {
      *interrupt_flag_ = true;
   }

   if( handle_interrupt_ != NULL )
   {
      (*handle_interrupt_)();
   }

   if( --abortcountdown_ == 0 )
   {
      fputs("Ipopt sighandler: Too many interrupt signals. Forcing termination.\n", stderr);
      exit(1);
   }
}

bool RegisterInterruptHandler(
   void        (*handle_interrupt)(void),
   bool*         interrupt_flag,
   unsigned int  abortlimit
)
{
   if( registered_handler )
   {
      return false;
   }
   registered_handler = true;
   abortcountdown_ = abortlimit;

   handle_interrupt_ = handle_interrupt;
   interrupt_flag_ = interrupt_flag;

#ifdef IPOPT_HAS_SIGACTION
   struct sigaction sa;
   sa.sa_handler = &sighandler;
   sa.sa_flags = SA_RESTART;
   sigfillset(&sa.sa_mask);
   if( sigaction(SIGINT, &sa, NULL) == -1 )
   {
      return false;
   }
   if( sigaction(SIGHUP, &sa, NULL) == -1 )
   {
      return false;
   }

#else
   signal(SIGINT, sighandler);
   signal(SIGTERM, sighandler);
   signal(SIGABRT, sighandler);

#endif

   return true;
}

bool UnregisterInterruptHandler(void)
{
   if( !registered_handler )
   {
      return false;
   }

#ifdef IPOPT_HAS_SIGACTION
   struct sigaction sa;
   sa.sa_handler = SIG_DFL;
   sa.sa_flags = SA_RESTART;
   sigfillset(&sa.sa_mask);
   if( sigaction(SIGINT, &sa, NULL) == -1 )
   {
      return false;
   }
   if( sigaction(SIGHUP, &sa, NULL) == -1 )
   {
      return false;
   }

#else
   signal(SIGINT, SIG_DFL);
   signal(SIGTERM, SIG_DFL);
   signal(SIGABRT, SIG_DFL);

#endif

   registered_handler = false;

   return true;
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
   va_list ap;
   va_start(ap, format);
   int ret;
   ret = vsnprintf(str, size, format, ap);
   va_end(ap);
   return ret;
}

} //namespace Ipopt
