// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpRestoFilterConvCheck.hpp"
#include "IpCompoundVector.hpp"
#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()
      :
      orig_filter_line_search_(NULL)
  {
    DBG_START_FUN("RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()",
		  dbg_verbosity);

  }

  RestoFilterConvergenceCheck::~RestoFilterConvergenceCheck()
  {
    DBG_START_FUN("RestoFilterConvergenceCheck::RestoFilterConvergenceCheck()",
		  dbg_verbosity);
  }

  void
  RestoFilterConvergenceCheck::SetOrigFilterLineSearch
  (const FilterLineSearch& orig_filter_line_search)
  {
    orig_filter_line_search_ = &orig_filter_line_search;
  }

  bool RestoFilterConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    DBG_ASSERT(orig_filter_line_search_ && "Need to call RestoFilterConvergenceCheck::SetOrigFilterLineSearch before Initialize");
    Number value = 0.0;

    // Check for the algorithm options
    if (options.GetNumericValue("kappa_resto", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0,
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"kappa_resto\": This value must be larger than 0 and less than 1.");
      kappa_resto_ = value;
    }
    else {
      kappa_resto_ = 0.9;
    }

    return OptimalityErrorConvergenceCheck::InitializeImpl(options, prefix);
  }

  ConvergenceCheck::ConvergenceStatus
  RestoFilterConvergenceCheck::CheckConvergence()
  {
    // First check if the point is now acceptable for the outer filter
    ConvergenceStatus status;

    // Get pointers to the Original NLP objects
    const RestoIpoptNLP* resto_ipopt_nlp =
      dynamic_cast<const RestoIpoptNLP*>(&IpNLP());
    DBG_ASSERT(resto_ipopt_nlp);

    SmartPtr<IpoptData> orig_ip_data = &resto_ipopt_nlp->OrigIpData();
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq =
      &resto_ipopt_nlp->OrigIpCq();

    // set the trial point for the original problem
    SmartPtr<const Vector> x = IpData().curr_x();
    const CompoundVector* cx =
      dynamic_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(cx);

    SmartPtr<const Vector> x_only = cx->GetComp(0);
    orig_ip_data->SetTrialPrimalVariablesFromPtr(x_only, IpData().curr_s());

    // Calculate the f and theta for the original problem
    Number orig_trial_theta = orig_ip_cq->trial_constraint_violation();
    Number orig_trial_barr = orig_ip_cq->trial_barrier_obj();
    Number orig_curr_theta = orig_ip_cq->curr_constraint_violation();

    // check acceptability to the filter
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_theta = %g, orig_trial_theta = %g, orig_trial_barr = %g\n",
                   orig_curr_theta, orig_trial_theta, orig_trial_barr);

    // ToDo: In the following we might want to be more careful with the lower bound
    Number orig_theta_max = Max(kappa_resto_*orig_curr_theta, orig_ip_data ->epsilon_tol());
    if (orig_trial_theta > orig_theta_max) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point does not provide sufficient reduction w.r.t the original theta.\n");
      status = CONTINUE;
    }
    else if (!orig_filter_line_search_->IsAcceptableToCurrentFilter(orig_trial_barr, orig_trial_theta)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point is not acceptable to the original filter.\n");
      status = CONTINUE;
    }
    else if (!orig_filter_line_search_->IsAcceptableToCurrentIterate(orig_trial_barr, orig_trial_theta) ) {
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

    // If the point is not yet acceptable to the filter, check if the problem
    // is maybe locally infeasible

    if (status==CONTINUE) {

      status = OptimalityErrorConvergenceCheck::CheckConvergence();
      if (status == CONVERGED) {
        THROW_EXCEPTION(LOCALLY_INFEASIBILE,
                        "Restoration phase converged to a point of local infeasibility");
      }
    }
    return status;
  }

} // namespace Ipopt
