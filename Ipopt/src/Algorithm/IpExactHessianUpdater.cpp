// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.cpp

#include "IpExactHessianUpdater.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  bool ExactHessianUpdater::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  void ExactHessianUpdater::UpdateHessian()
  {
    DBG_START_METH("ExactHessianUpdater::UpdateHessian",
                   dbg_verbosity);

    IpData().Set_W(IpCq().curr_exact_hessian());
  }


} // namespace Ipopt
