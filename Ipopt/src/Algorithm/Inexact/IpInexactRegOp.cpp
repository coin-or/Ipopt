// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-05

#include "IpInexactRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpInexactSearchDirCalc.hpp"
#include "IpInexactDoglegNormal.hpp"
#include "IpInexactNewtonNormal.hpp"
#include "IpInexactPDSolver.hpp"
#include "IpInexactLSAcceptor.hpp"

namespace Ipopt
{

  void RegisterOptions_Inexact(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Inexact Step Computation");
    InexactSearchDirCalculator::RegisterOptions(roptions);
    InexactDoglegNormalStep::RegisterOptions(roptions);
    InexactNewtonNormalStep::RegisterOptions(roptions);
    InexactPDSolver::RegisterOptions(roptions);
    InexactLSAcceptor::RegisterOptions(roptions);
  }

} // namespace Ipopt
