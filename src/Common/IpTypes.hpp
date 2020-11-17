// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTYPES_HPP__
#define __IPTYPES_HPP__

#include "IpoptConfig.h"

namespace Ipopt
{
/** Type of all numbers */
#ifdef IPOPT_SINGLE
  typedef float Number;
#else
  typedef double Number;
#endif
/** Type of all indices of vectors, matrices etc */
typedef int Index;
/** Type of default integer */
typedef int Int;

} // namespace Ipopt

/* Type of Fortran integer translated into C */
typedef IPOPT_FORTRAN_INTEGER_TYPE ipfint;

#endif
