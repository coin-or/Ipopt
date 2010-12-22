// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-06-24
//             based on IpRestoFilterConvCheck.cpp

#include "IpRestoPenaltyConvCheck.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  RestoPenaltyConvergenceCheck::RestoPenaltyConvergenceCheck()
      :
      orig_penalty_ls_acceptor_(NULL)
  {
    DBG_START_FUN("RestoPenaltyConvergenceCheck::RestoPenaltyConvergenceCheck()",
                  dbg_verbosity);
  }

  RestoPenaltyConvergenceCheck::~RestoPenaltyConvergenceCheck()
  {
    DBG_START_FUN("~RestoPenaltyConvergenceCheck::RestoPenaltyConvergenceCheck()",
                  dbg_verbosity);
  }

  void
  RestoPenaltyConvergenceCheck::SetOrigLSAcceptor
  (const BacktrackingLSAcceptor& orig_ls_acceptor)
  {
    orig_penalty_ls_acceptor_ = dynamic_cast<const PenaltyLSAcceptor*>(&orig_ls_acceptor);
    DBG_ASSERT(orig_penalty_ls_acceptor_);
  }

  void RestoPenaltyConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool RestoPenaltyConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    DBG_ASSERT(orig_penalty_ls_acceptor_ && "Need to call RestoPenaltyConvergenceCheck::SetOrigPenaltyLineSearch before Initialize");

    return RestoConvergenceCheck::InitializeImpl(options, prefix);
  }

  ConvergenceCheck::ConvergenceStatus
  RestoPenaltyConvergenceCheck::TestOrigProgress(Number orig_trial_barr,
      Number orig_trial_theta)
  {
    ConvergenceStatus status;

    if (!orig_penalty_ls_acceptor_->IsAcceptableToCurrentIterate(orig_trial_barr, orig_trial_theta, true) ) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point is not acceptable to the original current point.\n");
      status = CONTINUE;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Restoration found a point that provides sufficient reduction in"
                     " theta and is acceptable to the current penalty function.\n");
      status = CONVERGED;
    }

    return status;
  }

} // namespace Ipopt
