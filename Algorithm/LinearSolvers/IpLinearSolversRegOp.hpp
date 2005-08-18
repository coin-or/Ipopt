// Copyright (C) 2005, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLinearSolversRegOp.hpp 430 2005-08-10 00:19:54Z andreasw $
//
// Authors:  Carl Laird     IBM    2005-08-16

#ifndef __IPLINEARSOLVERSREGOP_HPP__
#define __IPLINEARSOLVERSREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_LinearSolvers(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
