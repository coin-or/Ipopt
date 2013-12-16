// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//
//           was originally IpRestoFilterConvCheck.hpp (rev 915)
//             separated by A Waechter IBM  2008-06-24

#include "IpRestoConvCheck.hpp"
#include "IpCompoundVector.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpRestoPhase.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  RestoConvergenceCheck::RestoConvergenceCheck()
  {
    DBG_START_FUN("RestoConvergenceCheck::RestoConvergenceCheck()",
                  dbg_verbosity);
  }

  RestoConvergenceCheck::~RestoConvergenceCheck()
  {
    DBG_START_FUN("~RestoConvergenceCheck::RestoConvergenceCheck()",
                  dbg_verbosity);
  }

  void RestoConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption(
      "required_infeasibility_reduction",
      "Required reduction of infeasibility before leaving restoration phase.",
      0.0, false, 1.0, true,
      0.9,
      "The restoration phase algorithm is performed, until a point is found "
      "that is acceptable to the filter and the infeasibility has been "
      "reduced by at least the fraction given by this option.");
    roptions->AddLowerBoundedIntegerOption(
      "max_resto_iter",
      "Maximum number of successive iterations in restoration phase.",
      0, 3000000,
      "The algorithm terminates with an error message if the number of "
      "iterations successively taken in the restoration phase exceeds this "
      "number.");
  }

  bool RestoConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("required_infeasibility_reduction", kappa_resto_, prefix);
    options.GetIntegerValue("max_iter", maximum_iters_, prefix);
    options.GetIntegerValue("max_resto_iter", maximum_resto_iters_, prefix);

    // The original constraint violation tolerance
    options.GetNumericValue("constr_viol_tol", orig_constr_viol_tol_, "");

    first_resto_iter_ = true;
    successive_resto_iter_ = 0;

    return OptimalityErrorConvergenceCheck::InitializeImpl(options, prefix);
  }

  ConvergenceCheck::ConvergenceStatus
  RestoConvergenceCheck::CheckConvergence(bool call_intermediate_callback /*= true*/)
  {
    // Get pointers to the Original NLP objects
    const RestoIpoptNLP* resto_ipopt_nlp =
      static_cast<const RestoIpoptNLP*>(&IpNLP());
    DBG_ASSERT(dynamic_cast<const RestoIpoptNLP*>(&IpNLP()));

    SmartPtr<IpoptData> orig_ip_data = &resto_ipopt_nlp->OrigIpData();
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq =
      &resto_ipopt_nlp->OrigIpCq();

    // set the trial point for the original problem
    SmartPtr<const Vector> x = IpData().curr()->x();
    const CompoundVector* cx =
      static_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(x)));
    SmartPtr<const Vector> s = IpData().curr()->s();
    const CompoundVector* cs =
      static_cast<const CompoundVector*>(GetRawPtr(s));
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(s)));
    DBG_ASSERT(cs->NComps() == 1);

    SmartPtr<IteratesVector> trial = orig_ip_data->curr()->MakeNewContainer();
    trial->Set_x(*cx->GetComp(0));
    trial->Set_s(*cs->GetComp(0));
    orig_ip_data->set_trial(trial);

    if (call_intermediate_callback) {
      // Check if user requested termination by calling the intermediate
      // user callback function
      AlgorithmMode mode = RestorationPhaseMode;
      // Gather the information also used in the iteration output
      Index iter = IpData().iter_count();
      Number inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
      Number inf_du = IpCq().curr_dual_infeasibility(NORM_MAX);
      Number mu = IpData().curr_mu();
      Number dnrm;
      if (IsValid(IpData().delta()) && IsValid(IpData().delta()->x()) && IsValid(IpData().delta()->s())) {
        dnrm = Max(IpData().delta()->x()->Amax(), IpData().delta()->s()->Amax());
      }
      else {
        // This is the first iteration - no search direction has been
        // computed yet.
        dnrm = 0.;
      }
      Number alpha_primal = IpData().info_alpha_primal();
      Number alpha_dual = IpData().info_alpha_dual();
      Number regu_x = IpData().info_regu_x();
      Number unscaled_f = orig_ip_cq->unscaled_trial_f();
      Index ls_count = IpData().info_ls_count();
      bool request_stop =
        !IpNLP().IntermediateCallBack(mode, iter, unscaled_f, inf_pr, inf_du,
                                      mu, dnrm, regu_x, alpha_dual,
                                      alpha_primal, ls_count,
                                      &IpData(), &IpCq());

      if (request_stop) {
        return ConvergenceCheck::USER_STOP;
      }
    }

    if (IpData().iter_count() >= maximum_iters_) {
      return ConvergenceCheck::MAXITER_EXCEEDED;
    }

    if (successive_resto_iter_ > maximum_resto_iters_) {
      Jnlst().Printf(J_WARNING, J_MAIN,
                     "More than %d successive iterations taken in restoration phase.\n",
                     maximum_resto_iters_);
      return ConvergenceCheck::MAXITER_EXCEEDED;
    }
    successive_resto_iter_++;

    // First check if the point is now acceptable for the outer filter
    ConvergenceStatus status;

    // Calculate the f and theta for the original problem
    Number orig_trial_theta = orig_ip_cq->trial_constraint_violation();
    Number orig_curr_theta = orig_ip_cq->curr_constraint_violation();

    // check acceptability to the filter
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_theta = %8.2e, orig_trial_theta = %8.2e\n",
                   orig_curr_theta, orig_trial_theta);

    // ToDo: In the following we might want to be more careful with the lower bound

    Number orig_curr_inf_pr = orig_ip_cq->curr_primal_infeasibility(NORM_MAX);
    Number orig_trial_inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_inf_pr = %8.2e, orig_trial_inf_pr = %8.2e\n",
                   orig_curr_inf_pr, orig_trial_inf_pr);


    Number orig_inf_pr_max = Max(kappa_resto_*orig_curr_inf_pr,
                                 Min(orig_ip_data->tol(),
                                     orig_constr_viol_tol_));
    if (kappa_resto_ == 0.) {
      orig_inf_pr_max = 0.;
    }

    if (first_resto_iter_) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "This is the first iteration - continue to take at least one step.\n");
      status = CONTINUE;
    }
    else if (orig_ip_cq->IsSquareProblem() &&
             orig_trial_inf_pr <=
             Min(orig_ip_data->tol(), orig_constr_viol_tol_)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Restoration phase found points satisfying feasibility tolerance in square problem.\n");
      status = CONVERGED;
    }
    else if (orig_trial_inf_pr > orig_inf_pr_max) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Point does not provide sufficient reduction w.r.t the original constraint violation (orig_inf_pr_max=%e).\n", orig_inf_pr_max);
      status = CONTINUE;
    }
    else {
      Number orig_trial_barr = orig_ip_cq->trial_barrier_obj();

      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "orig_trial_barr = %8.2e\n", orig_trial_barr);
      status = TestOrigProgress(orig_trial_barr, orig_trial_theta);
    }

    // If the point is not yet acceptable to the filter, check if the problem
    // is maybe locally infeasible

    if (status==CONTINUE) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Checking convergence for restoration phase problem...\n");
      status = OptimalityErrorConvergenceCheck::CheckConvergence(false);
      if (status == CONVERGED || status == CONVERGED_TO_ACCEPTABLE_POINT) {
        Number orig_trial_primal_inf =
          orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
        // ToDo make the factor in following line an option
        if (orig_trial_primal_inf <= 1e2*IpData().tol()) {
          //        if (orig_trial_primal_inf <= 1e2*orig_ip_data->tol()) {
          if (IpData().tol() > 1e-1*orig_ip_data->tol()) {
            // For once, we tighten the convergence tolerance for the
            // restoration phase problem in case the problem is only
            // very slightly infeasible.
            IpData().Set_tol(1e-2*IpData().tol());
            status = CONTINUE;
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "Tightening restoration phase tolerance to %e.\n",
                           IpData().tol());
            IpData().Append_info_string("!");
          }
          else {
            Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                           "Restoration phase converged to a feasible point that is\n"
                           "unacceptable to the filter for the original problem.\n");
            THROW_EXCEPTION(RESTORATION_CONVERGED_TO_FEASIBLE_POINT,
                            "Restoration phase converged to a feasible point that is "
                            "unacceptable to the filter for the original problem.");
          }
        }
        else {
          THROW_EXCEPTION(LOCALLY_INFEASIBLE,
                          "Restoration phase converged to a point of local infeasibility");
        }
      }
    }

    first_resto_iter_ = false;

    return status;
  }

} // namespace Ipopt
