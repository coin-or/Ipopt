// Copyright (C) 2007, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter         IBM        2007-06-01

#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpCGPenaltyRegOp.hpp"
#include "IpCGSearchDirCalc.hpp"
#include "IpCGPenaltyLSAcceptor.hpp"

namespace Ipopt
{

void RegisterOptions_CGPenalty(
   const SmartPtr<RegisteredOptions>& roptions
)
{
   roptions->SetRegisteringCategory("CG Penalty");
   CGSearchDirCalculator::RegisterOptions(roptions);
   CGPenaltyLSAcceptor::RegisterOptions(roptions);
   CGPenaltyCq::RegisterOptions(roptions);
}

} // namespace Ipopt
