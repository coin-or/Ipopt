// Copyright (C) 2004, International Business Machines and others.
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

  typedef double Number;
  typedef int Index;
  typedef int Int;

  enum ApplicationReturnStatus
  {
    Solve_Succeeded,
    Maximum_Iterations_Exceeded,
    Infeasible_Problem_Detected,
    Not_Enough_Degrees_Of_Freedom,
    Solve_Failed,
    Solved_To_Best_Possible_Precision,
    Solved_To_Acceptable_Level,
    NonIpopt_Exception_Thrown,
    Insufficient_Memory,
    Internal_Error
  };

} // namespace Ipopt

// external to the Ipopt namespace
typedef int ipfint;

#endif
