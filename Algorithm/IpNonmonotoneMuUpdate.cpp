// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpNonmonotoneMuUpdate.hpp"
#include "IpJournalist.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  NonmonotoneMuUpdate::NonmonotoneMuUpdate
  (const SmartPtr<LineSearch>& line_search,
   const SmartPtr<MuOracle>& mu_oracle)
      :
      MuUpdate(),
      linesearch_(line_search),
      mu_oracle_(mu_oracle)
  {
    DBG_ASSERT(IsValid(linesearch_));
    DBG_ASSERT(IsValid(mu_oracle_));
  }

  NonmonotoneMuUpdate::~NonmonotoneMuUpdate()
  {}

  bool NonmonotoneMuUpdate::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Number value;

    if (options.GetNumericValue("mu_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"mu_max\": This value must be larger than 0.");
      mu_max_ = value;
    }
    else {
      mu_max_ = 1e10;
    }

    if (options.GetNumericValue("mu_min", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < mu_max_, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"mu_min\": This value must be larger than 0 and less than mu_max.");
      mu_min_ = value;
    }
    else {
      mu_min_ = 0.1*IpData().epsilon_tol();
    }

    if (options.GetNumericValue("tau_min", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"tau_min\": This value must be between 0 and 1.");
      tau_min_ = value;
    }
    else {
      tau_min_ = 0.99;
    }

    if (options.GetNumericValue("tau_max", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"tau_max\": This value must be between 0 and 1.");
      tau_max_ = value;
    }
    else {
      tau_max_ = tau_min_;
    }

    if (options.GetNumericValue("nonmonotone_mu_refs_redfact", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0 && value < 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"nonmonotone_mu_refs_redfact\": This value must be between 0 and 1.");
      refs_red_fact_ = value;
    }
    else {
      refs_red_fact_ = 0.9999;
    }

    Index ivalue;
    if (options.GetIntegerValue("nonmonotone_mu_max_refs", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue >= 0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"nonmonotone_mu_max_refs\": This value must be non-negative.");
      num_refs_max_ = ivalue;
    }
    else {
      num_refs_max_ = 4;
    }

    if (options.GetIntegerValue("mu_never_fix", ivalue, prefix)) {
      mu_never_fix_ = (ivalue != 0);
    }
    else {
      mu_never_fix_ = false;
    }

    bool retvalue = mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                           options, prefix);

    refs_vals_.clear();
    check_if_no_bounds_ = false;
    no_bounds_ = false;
    fixed_mu_mode_ = false;

    // TODO do we need to initialize the linesearch object?

    return retvalue;
  }

  void NonmonotoneMuUpdate::UpdateBarrierParameter()
  {
    // of there are not bounds, we always return the minimum MU value
    if (!check_if_no_bounds_) {
      Index n_bounds = IpData().curr_z_L()->Dim() + IpData().curr_z_U()->Dim()
                       + IpData().curr_v_L()->Dim() + IpData().curr_v_U()->Dim();

      if( n_bounds==0 ) {
        no_bounds_ = true;
        IpData().Set_mu(mu_min_);
        IpData().Set_tau(tau_min_);
      }

      check_if_no_bounds_ = true;
    }

    if (no_bounds_)
      return;

    if (fixed_mu_mode_) {
      // if we are in the fixed mu mode, we need to check if the
      // current iterate is good enough to continue with the free mode
      bool sufficient_progress = CheckSufficientProgress();
      if (sufficient_progress) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Switching back to free mu mode.\n");
        fixed_mu_mode_ = false;
        RememberCurrentPointAsAccepted();
      }
      else {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Remaining in fixed mu mode.\n");
      }
    }
    else {
      bool sufficient_progress = CheckSufficientProgress();
      if (sufficient_progress) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Staying in free mu mode.\n");
        RememberCurrentPointAsAccepted();
      }
      else {
        fixed_mu_mode_ = true;

        // Set the new values for mu and tau and tell the linesearch
        // to reset its memory
        Number mu = NewFixedMu();
        Number tau = Compute_tau(mu);

        IpData().Set_mu(mu);
        IpData().Set_tau(tau);
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Switching to fixed mu mode with mu = %e and tau = %e.\n", mu, tau);
        linesearch_->Reset();
      }
    }

    if (!fixed_mu_mode_) {
      // Compute the new barrier parameter via the oracle
      Number mu = mu_oracle_->CalculateMu();

      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Barrier parameter mu computed by oracle is %e\n",
                     mu);

      // Apply safeguards if appropriate
      mu = Min(mu, mu_max_);
      mu = Max(mu, mu_min_);
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Barrier parameter mu after safeguards is %e\n",
                     mu);

      // Update the fraction-to-the-boundary rule parameter
      // TODO The first rule makes tau too small early on.
      //    Number tau = Max(tau_min_, 1.-mu);
      Number tau = Compute_tau(mu);
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Fraction-to-the-boundary parameter tau is %e\n",
                     tau);

      // Set the new values
      IpData().Set_mu(mu);
      IpData().Set_tau(tau);

      linesearch_->Reset();
    }
    else {
      IpData().Append_info_string("F");
    }
  }

  bool
  NonmonotoneMuUpdate::CheckSufficientProgress()
  {
    if (mu_never_fix_) return true;

    bool retval = true;

    Index num_refs = refs_vals_.size();
    if (num_refs >= num_refs_max_) {
      retval = false;
      Number curr_error = curr_norm_pd_system();
      std::list<Number>::iterator iter;
      for (iter = refs_vals_.begin(); iter != refs_vals_.end();
           iter++) {
        if ( curr_error <= refs_red_fact_*(*iter) ) {
          retval = true;
        }
      }
    }

    return retval;
  }

  void
  NonmonotoneMuUpdate::RememberCurrentPointAsAccepted()
  {
    Number curr_error = curr_norm_pd_system();
    Index num_refs = refs_vals_.size();
    if (num_refs >= num_refs_max_) {
      refs_vals_.pop_front();
    }
    refs_vals_.push_back(curr_error);

    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_BARRIER_UPDATE)) {
      Index num_refs = 0;
      std::list<Number>::iterator iter;
      for (iter = refs_vals_.begin(); iter != refs_vals_.end();
           iter++) {
        num_refs++;
        Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                       "pd system reference[%2d] = %.6e\n", num_refs, *iter);
      }
    }
  }

  Number
  NonmonotoneMuUpdate::Compute_tau(Number mu)
  {
    return Max(tau_min_, Min(1.-mu, tau_max_));
    //return tau_min_;
  }

  Number
  NonmonotoneMuUpdate::NewFixedMu()
  {
    DBG_ASSERT(refs_vals_.size()>0);
    std::list<Number>::iterator iter = refs_vals_.begin();
    Number min_ref = *iter;
    iter++;
    while (iter != refs_vals_.end()) {
      min_ref = Min(min_ref, *iter);
      iter++;
    }

    Number avrg_compl = IpCq().curr_avrg_compl();
    return Min(avrg_compl, 0.1 * min_ref);
  }

  Number
  NonmonotoneMuUpdate::curr_norm_pd_system()
  {
    Number dual_inf =
      IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_1);
    Number primal_inf =
      IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_1);
    Number complty =
      IpCq().curr_complementarity(0., IpoptCalculatedQuantities::NORM_1);

    // scale those values (to get the average)
    Index n_dual = IpData().curr_x()->Dim() + IpData().curr_s()->Dim();
    dual_inf /= (Number)n_dual;
    Index n_pri = IpData().curr_y_c()->Dim() + IpData().curr_y_d()->Dim();
    DBG_ASSERT(n_pri>0 || primal_inf==0.);
    if (n_pri>0) {
      primal_inf /= (Number)n_pri;
    }
    Index n_comp = IpData().curr_z_L()->Dim() + IpData().curr_z_U()->Dim() +
                   IpData().curr_v_L()->Dim() + IpData().curr_v_U()->Dim();
    DBG_ASSERT(n_comp>0 || complty==0.);
    if (n_comp>0) {
      complty /= (Number)n_comp;
    }

    Number norm_pd_system = primal_inf + dual_inf + complty;

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "In barrier update check:\n"
                   "  average primal infeasibility: %15.6e\n"
                   "    average dual infeasibility: %15.6e\n"
                   "       average complementarity: %15.6e\n"
                   "   scaled norm of pd equations: %15.6e\n",
                   primal_inf, dual_inf, complty, norm_pd_system);

    return norm_pd_system;
  }

} // namespace Ipopt
