// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTYPES_HPP__
#define __IPTYPES_HPP__

#include "config.h"

namespace Ipopt
{
  /** Type of all numbers */
  typedef double Number;
  /** Type of all indices of vectors, matrices etc */
  typedef int Index;
  /** Type of default integer */
  typedef int Int;

} // namespace Ipopt

// external to the Ipopt namespace
/** Type of Fortran integer translated into C */
typedef int ipfint;

#endif
