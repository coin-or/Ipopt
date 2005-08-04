// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpException.hpp 277 2005-05-04 20:54:12Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPINTERFACETYPES_HPP__
#define __IPINTERFACETYPES_HPP__

namespace Ipopt
{
  enum ApplicationReturnStatus
  {
    Solve_Succeeded,
    Maximum_Iterations_Exceeded,
    Solved_To_Best_Possible_Precision,
    Solved_To_Acceptable_Level,
    Infeasible_Problem_Detected,
    Restoration_Failed,
    Not_Enough_Degrees_Of_Freedom,
    Unrecoverable_Exception,
    NonIpopt_Exception_Thrown,
    Insufficient_Memory,
    Internal_Error
  };

};

#endif
