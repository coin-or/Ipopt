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

  DefineIpoptType(RestoFilterConvergenceCheck);

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

  void RestoFilterConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption("kappa_resto", "???", 
				     0.0, true, 1.0, true, 0.9);
  }

  bool RestoFilterConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    DBG_ASSERT(orig_filter_line_search_ && "Need to call RestoFilterConvergenceCheck::SetOrigFilterLineSearch before Initialize");
    options.GetNumericValue("kappa_resto", kappa_resto_, prefix);

    first_resto_iter_ = true;

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
    SmartPtr<const Vector> x = IpData().curr()->x();
    const CompoundVector* cx =
      dynamic_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(cx);

    SmartPtr<IteratesVector> trial = orig_ip_data->curr()->MakeNewContainer();
    trial->Set_x(*cx->GetComp(0));
    trial->Set_s(*IpData().curr()->s());
    orig_ip_data->set_trial(trial);

    // Calculate the f and theta for the original problem
    Number orig_trial_theta = orig_ip_cq->trial_constraint_violation();
    Number orig_trial_barr = orig_ip_cq->trial_barrier_obj();
    Number orig_curr_theta = orig_ip_cq->curr_constraint_violation();

    // check acceptability to the filter
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_theta = %8.2e, orig_trial_theta = %8.2e, orig_trial_barr = %8.2e\n",
                   orig_curr_theta, orig_trial_theta, orig_trial_barr);

    // ToDo: In the following we might want to be more careful with the lower bound
    Number orig_theta_max = Max(kappa_resto_*orig_curr_theta,
                                1.e1*Min(orig_ip_data->tol(),
                                         orig_ip_data->primal_inf_tol()));
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
    else if (!orig_filter_line_search_->IsAcceptableToCurrentIterate(orig_trial_barr, orig_trial_theta, true) ) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point is not acceptable to the original current point.\n");
      status = CONTINUE;
    }
    else if (first_resto_iter_) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "This is the first iteration - continue to take at least one step.\n");
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
        Number orig_trial_primal_inf =
          orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
        // ToDo make the factor in following line an option
        if (orig_trial_primal_inf <= 1e2*orig_ip_data->tol()) {
          THROW_EXCEPTION(RESTORATION_FAILED,
                          "Restoration phase converged to a point with small primal infeasibility");
        }
        else {
          THROW_EXCEPTION(LOCALLY_INFEASIBILE,
                          "Restoration phase converged to a point of local infeasibility");
        }
      }
    }

    first_resto_iter_ = false;

    return status;
  }

} // namespace Ipopt
