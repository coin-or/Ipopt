// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#ifndef __IPALGBUILDER_HPP__
#define __IPALGBUILDER_HPP__

#include "IpUtils.hpp"
#include "IpIpoptAlg.hpp"

namespace Ipopt
{
  /** Builder to create a complete IpoptAlg object.  This object
   *  contains all subelements (such as line search objects etc).  How
   *  the resulting IpoptAlg object is build can be influenced by the
   *  options. */
  SmartPtr<IpoptAlgorithm> AlgorithmBuilder(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix);

} // namespace Ipopt

#endif
