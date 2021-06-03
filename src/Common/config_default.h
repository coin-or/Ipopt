
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_ipopt_default.h */
#ifndef IPOPTLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
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
/* #define IPOPT_CHECKLEVEL 0 */

/* Define to the debug verbosity level (0 is no output) */
/* #define IPOPT_VERBOSITY 0 */

/* If defined, the Ampl Solver Library is available. */
/* #undef IPOPT_HAS_ASL 1 */

/* If defined, the LAPACK Library is available. */
#define IPOPT_HAS_LAPACK 1

/* If defined, the HSL library is available. */
/* #undef IPOPT_HAS_HSL 1 */

/* If defined, the MUMPS library is available. */
/* #undef IPOPT_HAS_MUMPS */

/* Define to 1 if the linear solver loader should be compiled to allow dynamic
   loading of shared libraries with linear solvers */
/* #undef IPOPT_HAS_LINEARSOLVERLOADER */

/* Define to 1 if you are using Pardiso from MKL */
/* #undef IPOPT_HAS_PARDISO_MKL */

/* Define to 1 if SPRAL is available */
/* #undef IPOPT_HAS_SPRAL */

/* Define to 1 if WSMP is available */
/* #undef IPOPT_HAS_WSMP */

/* Define to be the name of C-function for Inf check */
#ifdef _MSC_VER
#define IPOPT_C_FINITE _finite
#else
#define IPOPT_C_FINITE std::isfinite
#endif

#define IPOPT_BLAS_FUNC(name,NAME)    F77_FUNC(name,NAME)
#define IPOPT_LAPACK_FUNC(name,NAME)  F77_FUNC(name,NAME)
#define IPOPT_PARDISO_FUNC(name,NAME) F77_FUNC(name,NAME)
#define IPOPT_WSMP_FUNC(name,NAME)    F77_FUNC(name,NAME)
#define IPOPT_WSMP_FUNC_(name,NAME)   F77_FUNC_(name,NAME)
