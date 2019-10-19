
/***************************************************************************/
/*           HERE DEFINE THE PROJECT SPECIFIC PUBLIC MACROS                */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Version number of project */
#define IPOPT_VERSION      "3.13.0"

/* Major Version number of project */
#define IPOPT_VERSION_MAJOR      3

/* Minor Version number of project */
#define IPOPT_VERSION_MINOR     13

/* Release Version number of project */
#define IPOPT_VERSION_RELEASE    0

/* Define to the C type corresponding to Fortran INTEGER */
#ifndef FORTRAN_INTEGER_TYPE
#define FORTRAN_INTEGER_TYPE int
#endif

#ifndef IPOPTLIB_EXPORT
#ifdef _WIN32
/* assuming we link against an Ipopt DLL */
#define IPOPTLIB_EXPORT __declspec(dllimport)
#else
#define IPOPTLIB_EXPORT
#endif
