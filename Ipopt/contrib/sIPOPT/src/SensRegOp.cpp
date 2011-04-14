// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
// modified for sIPOPT by Hans Pirnay, 2009-07-22

#include "IpRegOptions.hpp"
#include "SensApplication.hpp"

namespace Ipopt
{
  void RegisterOptions_sIPOPT(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Uncategorized");
    SensApplication::RegisterOptions(roptions);
  }

} // namespace Ipopt
