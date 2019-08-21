// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-14


#ifndef __SENSCUTILS_HPP__
#define __SENSCUTILS_HPP__

#include "IpUtils.hpp"
#include "IpAlgTypes.hpp"
#include "IpReturnCodes.hpp"

#include <string>

namespace Ipopt
{

  /** This header file provides some definitions used throughout the program. */


  enum SensAlgorithmExitStatus{
    SOLVE_SUCCESS,
    FATAL_ERROR
  };

  SIPOPTLIB_EXPORT Index AsIndexMax(Index length, const Index* x, Index Incr);

  SIPOPTLIB_EXPORT Index AsIndexSum(Index length, const Index* x, Index Incr);

  SIPOPTLIB_EXPORT void append_Index(std::string& str, Index idx);

  SIPOPTLIB_EXPORT SolverReturn AppReturn2SolverReturn(ApplicationReturnStatus ipopt_retval);
}

#endif
