// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMonotoneMuUpdate.hpp"
#include "IpJournalist.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  DefineIpoptType(MonotoneMuUpdate);

  MonotoneMuUpdate::MonotoneMuUpdate(const SmartPtr<LineSearch>& linesearch)
      :
      MuUpdate(),
      linesearch_(linesearch),
      initialized_(false)
  {
    DBG_START_METH("MonotoneMuUpdate::MonotoneMuUpdate",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(linesearch_));
  }

  MonotoneMuUpdate::~MonotoneMuUpdate()
  {
    DBG_START_METH("MonotoneMuUpdate::~MonotoneMuUpdate",
                   dbg_verbosity);
  }

  void MonotoneMuUpdate::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption("mu0", "initial value for the barrier parameter, mu",
                                          0.0, true, 0.1);
    //     roptions->AddLowerBoundedNumberOption("kappa_epsilon", "???",
    // 					  0.0, true, 10.0);
    //     roptions->AddBoundedNumberOption("kappa_mu", "???",
    // 				     0.0, true, 1.0, true, 0.2);
    //     roptions->AddBoundedNumberOption("theta_mu", "???",
    // 				     1.0, true, 2.0, true, 1.5);
    //     roptions->AddBoundedNumberOption("tau_min", "???",
    //  				     0.0, true, 1.0, true, 0.99);
  }

  bool MonotoneMuUpdate::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("mu0", mu0_, prefix);
    options.GetNumericValue("kappa_epsilon", kappa_epsilon_, prefix);
    options.GetNumericValue("kappa_mu", kappa_mu_, prefix);
    options.GetNumericValue("theta_mu", theta_mu_, prefix);
    options.GetNumericValue("tau_min", tau_min_, prefix);

    IpData().Set_mu(mu0_);
    Number tau = Max(tau_min_, 1.0 - mu0_);
    IpData().Set_tau(tau);

    initialized_ = false;

    //TODO we need to clean up the mu-update for the restoration phase
    if (prefix=="resto.") {
      first_iter_resto_ = true;
    }
    else {
      first_iter_resto_ = false;
    }

    return true;
  }

  void MonotoneMuUpdate::UpdateBarrierParameter()
  {
    Number mu = IpData().curr_mu();
    Number tau = IpData().curr_tau();

    Number sub_problem_error = IpCq().curr_barrier_error();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Optimality Error for Barrier Sub-problem = %e\n",
                   sub_problem_error);
    Number kappa_eps_mu = kappa_epsilon_ * mu;

    bool done = false;
    bool tiny_step_flag = IpData().tiny_step_flag();
    while ((sub_problem_error <= kappa_eps_mu || tiny_step_flag)
           && !done && !first_iter_resto_) {
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "  sub_problem_error < kappa_eps * mu (%e)\n", kappa_eps_mu);

      // Compute the new values for mu and tau
      Number new_mu;
      Number new_tau;
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Updating mu=%e and tau=%e to ", mu, tau);
      CalcNewMuAndTau(new_mu, new_tau);
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "new_mu=%e and new_tau=%e\n", new_mu, new_tau);
      bool mu_changed = (mu != new_mu);
      if (!mu_changed && tiny_step_flag) {
        THROW_EXCEPTION(TINY_STEP_DETECTED,
                        "Problem solved to best possible numerical accuracy");
      }

      // Set the new values for mu and tau
      IpData().Set_mu(new_mu);
      IpData().Set_tau(new_tau);
      mu = new_mu;
      tau = new_tau;

      // If this is the first iteration, we want to check if we can
      // decrease mu even more
      if (initialized_) {
        done = true;
      }
      else {
        sub_problem_error = IpCq().curr_barrier_error();
        kappa_eps_mu = kappa_epsilon_ * mu;
        done = (sub_problem_error > kappa_eps_mu);
      }

      // Reset the line search
      if (done && mu_changed) {
        linesearch_->Reset();
      }

      tiny_step_flag = false;
    }

    first_iter_resto_ = false;
    initialized_ = true;
  }

  void MonotoneMuUpdate::CalcNewMuAndTau(Number &new_mu,
                                         Number &new_tau)
  {
    // update the barrier parameter
    Number mu = IpData().curr_mu();
    Number tol = IpData().tol();
    Number compl_inf_tol = IpData().compl_inf_tol();

    new_mu = Min( kappa_mu_*mu, pow(mu, theta_mu_) );
    new_mu = Max(new_mu, Min(tol, compl_inf_tol)/10.);

    // update the fraction to the boundary parameter
    new_tau = Max(tau_min_, 1.-new_mu);
  }

} // namespace Ipopt
