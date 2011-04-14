// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//

#ifndef __SENSREGOP_HPP__
#define __SENSREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_sIPOPT(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
