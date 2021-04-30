// Copyright (C) 2020 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef __IPTYPES_H__
#define __IPTYPES_H__

#include "IpoptConfig.h"

/// macro to declare symbols as deprecated
/// @since 3.14.0
#ifndef IPOPT_DEPRECATED
#if defined(_MSC_VER)
#  define IPOPT_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__)
#  define IPOPT_DEPRECATED __attribute__ ((deprecated))
#else
#  define IPOPT_DEPRECATED
#endif
#endif

/* Type of Fortran integer translated into C */
typedef IPOPT_FORTRAN_INTEGER_TYPE ipfint;

/** Type for floating-point numbers
 * @since 3.14.0
 */
#ifdef IPOPT_SINGLE
typedef float ipnumber;
#else
typedef double ipnumber;
#endif

/** Type of all indices of vectors, matrices etc
 * @since 3.14.0
 */
typedef int ipindex;

#endif
