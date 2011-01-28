// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#ifndef __IPLINEARSOLVERSREGOP_HPP__
#define __IPLINEARSOLVERSREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
