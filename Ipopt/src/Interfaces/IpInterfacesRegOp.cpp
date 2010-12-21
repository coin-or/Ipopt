// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
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
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptApplication::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    TNLPAdapter::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
