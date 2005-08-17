// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMa27TSolverInterface.hpp 430 2005-08-10 00:19:54Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpLinearSolversRegOp.hpp"
#include "IpRegOptions.hpp"
#include "IpMa27TSolverInterface.hpp"

namespace Ipopt
{

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions)
  {
    Ma27TSolverInterface::RegisterOptions(roptions);
  }

} // namespace Ipopt
