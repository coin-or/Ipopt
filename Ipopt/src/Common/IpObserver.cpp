// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpObserver.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG_OBSERVER
  const Index Observer::dbg_verbosity = 0;
  const Index Subject::dbg_verbosity = 0;
#endif
} // namespace Ipopt
