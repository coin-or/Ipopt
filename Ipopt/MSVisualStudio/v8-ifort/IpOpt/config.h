#include "config_default.h"

/* If defined, the METIS library is available. */
#define COIN_HAS_METIS 1

/* If defined, the MUMPS Library is available. */
#define COIN_HAS_MUMPS 1

#undef COIN_HAS_HSL

/* Define to 1 if the linear solver loader should be compiled to allow dynamic
   loading of shared libaries with linear solvers */
#define HAVE_LINEARSOLVERLOADER 1

/* Define to 1 if Pardiso is available */
#undef HAVE_PARDISO

#define HAVE_CSTDDEF 1

#define SHAREDLIBEXT "dll"

#define __str__(s) #s
#define __xstr__(s) __str__(s)

#define IPOPT_RESOURCE_VERSION IPOPT_VERSION_MAJOR,IPOPT_VERSION_MINOR,IPOPT_VERSION_RELEASE,0
