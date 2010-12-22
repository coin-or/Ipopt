// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//
//           A Waechter: moved most code to IpRestoConvCheck.cpp 2008-06-24

#include "IpRestoFilterConvCheck.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()
      :
      orig_filter_ls_acceptor_(NULL)
  {
    DBG_START_FUN("RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()",
                  dbg_verbosity);
  }

  RestoFilterConvergenceCheck::~RestoFilterConvergenceCheck()
  {
    DBG_START_FUN("~RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()",
                  dbg_verbosity);
  }

  void
  RestoFilterConvergenceCheck::SetOrigLSAcceptor
  (const BacktrackingLSAcceptor& orig_ls_acceptor)
  {
    orig_filter_ls_acceptor_ = dynamic_cast<const FilterLSAcceptor*>(&orig_ls_acceptor);
    DBG_ASSERT(orig_filter_ls_acceptor_);
  }

  void RestoFilterConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool RestoFilterConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    DBG_ASSERT(orig_filter_ls_acceptor_ && "Need to call RestoFilterConvergenceCheck::SetOrigFilterLineSearch before Initialize");

    return RestoConvergenceCheck::InitializeImpl(options, prefix);
  }

  ConvergenceCheck::ConvergenceStatus
  RestoFilterConvergenceCheck::TestOrigProgress(Number orig_trial_barr,
      Number orig_trial_theta)
  {
    ConvergenceStatus status;

    if (!orig_filter_ls_acceptor_->IsAcceptableToCurrentFilter(orig_trial_barr, orig_trial_theta)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point is not acceptable to the original filter.\n");
      status = CONTINUE;
    }
    else if (!orig_filter_ls_acceptor_->IsAcceptableToCurrentIterate(orig_trial_barr, orig_trial_theta, true) ) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point is not acceptable to the original current point.\n");
      status = CONTINUE;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Restoration found a point that provides sufficient reduction in"
                     " theta and is acceptable to the current filter.\n");
      status = CONVERGED;
    }

    return status;
  }

} // namespace Ipopt
