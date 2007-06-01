// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAlgorithmRegOp.cpp 984 2007-05-30 23:03:22Z andreasw $
//
// Authors:  Andreas Waechter         IBM        2007-06-01

#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpCGSearchDirCalc.hpp"
#include "IpCGPenaltyLSAcceptor.hpp"

namespace Ipopt
{

  void RegisterOptions_CGPenalty(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Line Search");
    CGSearchDirCalculator::RegisterOptions(roptions);
    CGPenaltyLSAcceptor::RegisterOptions(roptions);
  }

} // namespace Ipopt
