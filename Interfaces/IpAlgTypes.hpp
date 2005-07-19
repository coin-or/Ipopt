// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpException.hpp 277 2005-05-04 20:54:12Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPALGTYPES_HPP__
#define __IPALGTYPES_HPP__

#include "IpTypes.hpp"
#include "IpException.hpp"

namespace Ipopt
{

  /**@name Enumerations */
  //@{
  /** enum for the return from the optimize algorithm
   *  (obviously we need to add more) */
  enum SolverReturn {
    SUCCESS,
    MAXITER_EXCEEDED,
    STOP_AT_TINY_STEP,
    STOP_AT_ACCEPTABLE_POINT,
    LOCAL_INFEASIBILITY,
    RESTORATION_FAILURE
  };
  //@}

  /** Some exceptions used in multiple places */
  DECLARE_STD_EXCEPTION(LOCALLY_INFEASIBLE);
  DECLARE_STD_EXCEPTION(TOO_FEW_DOF);
  DECLARE_STD_EXCEPTION(TINY_STEP_DETECTED);
  DECLARE_STD_EXCEPTION(ACCEPTABLE_POINT_REACHED);

};

#endif
