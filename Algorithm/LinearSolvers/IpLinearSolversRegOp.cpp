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
#ifdef HAVE_MA57
# include "IpMa57TSolverInterface.hpp"
#endif
#ifdef HAVE_PARDISO
# include "IpPardisoSolverInterface.hpp"
#endif

namespace Ipopt
{

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions)
  {
#ifdef HAVE_MA27

    roptions->SetRegisteringCategory("MA27 Linear Solver");
    Ma27TSolverInterface::RegisterOptions(roptions);
#endif

#ifdef HAVE_MA57

    roptions->SetRegisteringCategory("MA57 Linear Solver");
    Ma57TSolverInterface::RegisterOptions(roptions);
#endif

#ifdef HAVE_PARDISO

    PardisoSolverInterface::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Pardiso Linear Solver");
#endif

    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
