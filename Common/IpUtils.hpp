// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPUTILS_HPP__
#define __IPUTILS_HPP__

// Standard Ip Include Files
#include "IpTypes.hpp"
#include "IpDebug.hpp"

#ifndef MY_C_FINITE
# define FiniteNumber finite
#else
# define FiniteNumber MY_C_FINITE
#endif

namespace Ipopt
{

  inline ipfint Max(ipfint a, ipfint b)
  {
    return ((a) > (b) ? (a) : (b));
  }

  inline Number Max(Number a, Number b)
  {
    return ((a) > (b) ? (a) : (b));
  }

  inline Number Max(Number a, Number b, Number c)
  {
    Number max = Max(a,b);
    max = Max(max, c);
    return max;
  }

  inline Number Max(Number a, Number b, Number c, Number d)
  {
    Number max = Max(a, b, c);
    max = Max(max, d);
    return max;
  }

  inline Number Min(Number a, Number b)
  {
    return ((a) < (b) ? (a) : (b));
  }

  inline Number Min(Number a, Number b, Number c)
  {
    Number min = Min(a,b);
    min = Min(min, c);
    return min;
  }

  inline Number Min(Number a, Number b, Number c, Number d)
  {
    Number min = Min(a, b, c);
    min = Min(min, d);
    return min;
  }

} //namespace Ipopt

#endif
