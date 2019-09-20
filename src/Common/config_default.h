
/* include the COIN-OR-wide system specific configure header */
/* this file lives in BuildTools/headers, which may not be around
 * however, with the normal autotools-based build system, we should
 * never get into a situation where this code should be compiled
 */
#warning "Ipopt-internal config_default.h is compiled. This is likely going to fail."
/* #include "configall_system.h" */

/* this needs to come before the include of config_ipopt_default.h */
#ifndef IPOPTLIB_EXPORT
#ifdef _WIN32
/* assuming we build an Ipopt DLL */
#define IPOPTLIB_EXPORT __declspec(dllexport)
#else
#define IPOPTLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_ipopt_default.h"

/***************************************************************************/
/*        HERE DEFINE THE PROJECT SPECIFIC PRIVATE MACROS                  */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to the debug sanity check level (0 is no test) */
#define COIN_IPOPT_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_IPOPT_VERBOSITY 0

/* If defined, the Ampl Solver Library is available. */
#define COIN_HAS_ASL 1

/* If defined, the LAPACK Library is available. */
#define COIN_HAS_LAPACK 1

/* If defined, the HSL library is available. */
#define COIN_HAS_HSL 1

/* If defined, the MUMPS library is available. */
/* #undef COIN_HAS_MUMPS */

/* Define to 1 if the linear solver loader should be compiled to allow dynamic
   loading of shared libraries with linear solvers */
/* #undef HAVE_LINEARSOLVERLOADER */

/* Define to 1 if Pardiso is available */
/* #undef HAVE_PARDISO */

/* Define to 1 if you are using Pardiso from MKL */
/* #undef HAVE_PARDISO_MKL */

/* Define to 1 if you are using the parallel version of Pardiso */
/* #undef HAVE_PARDISO_PARALLEL */

/* Define to 1 if WSMP is available */
/* #undef HAVE_WSMP */

/* Define to the C type corresponding to Fortran INTEGER */
#ifndef FORTRAN_INTEGER_TYPE
#define FORTRAN_INTEGER_TYPE int
#endif

#ifdef _MSC_VER
/* Define to be the name of C-function for Inf check */
#define COIN_C_FINITE _finite
#endif
