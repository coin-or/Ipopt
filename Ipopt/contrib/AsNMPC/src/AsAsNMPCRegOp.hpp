// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//

#ifndef __ASASNMPCREGOP_HPP__
#define __ASASNMPCREGOP_HPP__

#include "IpSmartPtr.hpp"

namespace Ipopt
{
  class RegisteredOptions;

  void RegisterOptions_AsNMPC(const SmartPtr<RegisteredOptions>& roptions);

} // namespace Ipopt

#endif
