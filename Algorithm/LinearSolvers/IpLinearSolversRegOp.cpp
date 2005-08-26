// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpLinearSolversRegOp.hpp"
#include "IpRegOptions.hpp"

#ifdef HAVE_MA27
# include "IpMa27TSolverInterface.hpp"
#endif

namespace Ipopt
{

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("MA27 Linear Solver");
#ifdef HAVE_MA27

    Ma27TSolverInterface::RegisterOptions(roptions);
#endif

    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
