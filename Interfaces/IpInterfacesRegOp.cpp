// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMa27TSolverInterface.hpp 430 2005-08-10 00:19:54Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpInterfacesRegOp.hpp"
#include "IpRegOptions.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLPAdapter.hpp"

namespace Ipopt
{

  void RegisterOptions_Interfaces(const SmartPtr<RegisteredOptions>& roptions)
  {
    IpoptApplication::RegisterOptions(roptions);
    TNLPAdapter::RegisterOptions(roptions);
  }

} // namespace Ipopt
