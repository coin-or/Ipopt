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

  bool MonotoneMuUpdate::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    Number value= 0.0;

    if (options.GetNumericValue("mu0", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"mu0\": This value must be larger than 0.");
      mu0_ = value;
    }
    else {
      mu0_ = 0.1;
    }

    if (options.GetNumericValue("kappa_epsilon", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"kappa_epsilon\": This value must be larger than 0.");
      kappa_epsilon_ = value;
    }
    else {
      kappa_epsilon_ = 10.0;
    }

    if (options.GetNumericValue("kappa_mu", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"kappa_mu\": This value must be between 0 and 1.");
      kappa_mu_ = value;
    }
    else {
      kappa_mu_ = 0.2;
    }

    if (options.GetNumericValue("theta_mu", value, prefix)) {
      ASSERT_EXCEPTION(value > 1.0 && value < 2.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"theta_mu\": This value must be between 1 and 2.");
      theta_mu_ = value;
    }
    else {
      theta_mu_ = 1.5;
    }

    if (options.GetNumericValue("tau_min", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"tau_min\": This value must be between 0 and 1.");
      tau_min_ = value;
    }
    else {
      tau_min_ = 0.99;
    }

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
    while (sub_problem_error <= kappa_eps_mu && !done && !first_iter_resto_) {
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
    }

    first_iter_resto_ = false;
    initialized_ = true;
  }

  void MonotoneMuUpdate::CalcNewMuAndTau(Number &new_mu,
                                         Number &new_tau)
  {
    // update the barrier parameter
    Number mu = IpData().curr_mu();
    Number eps_tol = IpData().epsilon_tol();

    new_mu = Min( kappa_mu_*mu, pow(mu, theta_mu_) );
    new_mu = Max(new_mu, eps_tol/10);

    // update the fraction to the boundary parameter
    new_tau = Max(tau_min_, 1.-new_mu);
  }

} // namespace Ipopt
