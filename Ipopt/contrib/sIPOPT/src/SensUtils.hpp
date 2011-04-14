// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-14


#ifndef __SENSCUTILS_HPP__
#define __SENSCUTILS_HPP__

#include "IpUtils.hpp"
#include <string>
#include "IpAlgTypes.hpp"
#include "IpReturnCodes.hpp"

namespace Ipopt
{

  /** This header file provides some definitions used throughout the program. */


  enum SensAlgorithmExitStatus{
    SOLVE_SUCCESS,
    FATAL_ERROR
  };

  Index AsIndexMax(Index length, const Index* x, Index Incr);

  Index AsIndexSum(Index length, const Index* x, Index Incr);

  void append_Index(std::string& str, Index idx);

  SolverReturn AppReturn2SolverReturn(ApplicationReturnStatus ipopt_retval);
}

#endif
