// Copyright (C) 2020 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef __IPTYPES_H__
#define __IPTYPES_H__

#include "IpoptConfig.h"

/* Type of Fortran integer translated into C */
typedef IPOPT_FORTRAN_INTEGER_TYPE ipfint;

/** Type for floating-point numbers
 *
 * Must be the same as Ipopt::Number.
 */
#ifdef IPOPT_SINGLE
typedef float ipnumber;
#else
typedef double ipnumber;
#endif

#endif
