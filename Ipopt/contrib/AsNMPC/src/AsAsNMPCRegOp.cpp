// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
// modified for AsNMPC by Hans Pirnay, 2009-07-22

#include "IpRegOptions.hpp"
#include "AsNMPCApplication.hpp"

namespace Ipopt
{
  void RegisterOptions_AsNMPC(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Uncategorized");
    NmpcApplication::RegisterOptions(roptions);
  }

} // namespace Ipopt
