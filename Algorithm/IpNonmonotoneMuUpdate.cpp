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
   const SmartPtr<MuOracle>& free_mu_oracle,
   const SmartPtr<MuOracle>& fix_mu_oracle)
      :
      MuUpdate(),
      linesearch_(line_search),
      free_mu_oracle_(free_mu_oracle),
      fix_mu_oracle_(fix_mu_oracle),
      filter_(2)
  {
    DBG_ASSERT(IsValid(linesearch_));
    DBG_ASSERT(IsValid(free_mu_oracle_));
    // fix_mu_oracle may be NULL
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
      ASSERT_EXCEPTION(value > 0.0 && value <= 1.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"tau_max\": This value must be between 0 and 1.");
      tau_max_ = value;
    }
    else {
      tau_max_ = tau_min_;
    }

    if (options.GetNumericValue("mu_safeguard_exp", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"mu_safeguard_exp\": This value must be non-negative.");
      mu_safeguard_exp_ = value;
    }
    else {
      mu_safeguard_exp_ = 0.;
    }

    if (options.GetNumericValue("mu_safeguard_factor", value, prefix)) {
      ASSERT_EXCEPTION(value >= 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"mu_safeguard_factor\": This value must be non-negative.");
      mu_safeguard_factor_ = value;
    }
    else {
      mu_safeguard_factor_ = 0.;
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

    if (options.GetIntegerValue("adaptive_globalization", ivalue, prefix)) {
      adaptive_globalization_ = ivalue;
    }
    else {
      adaptive_globalization_ = 1;
    }

    if (options.GetIntegerValue("restore_accepted_iterate", ivalue, prefix)) {
      restore_accepted_iterate_ = (ivalue != 0);
    }
    else {
      restore_accepted_iterate_ = false;
    }

    bool retvalue = free_mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(),
                    IpCq(), options, prefix);
    if (!retvalue) {
      return retvalue;
    }

    if (IsValid(fix_mu_oracle_)) {
      retvalue = fix_mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(),
                                            IpCq(), options, prefix);
      if (!retvalue) {
        return retvalue;
      }
    }

    if (options.GetNumericValue("fixed_mu_avrg_factor", value, prefix)) {
      ASSERT_EXCEPTION(value > 0.0, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"fixed_mu_avrg_factor\": This value must be larger than 0.");
      fixed_mu_avrg_factor_ = value;
    }
    else {
      fixed_mu_avrg_factor_ = 0.8;
    }

    // ToDo combine the following with MonotoneMuUpdate
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

    if (options.GetIntegerValue("nonmonotone_kkt_norm", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=1 && ivalue<=3, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"nonmonotone_kkt_norm\": This value must be between 1 and 3.");
      nonmonotone_kkt_norm_ = ivalue;
    }
    else if (options.GetIntegerValue("quality_function_norm", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=1 && ivalue<=4, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"quality_function_norm\": This value must be between 1 and 3.");
      nonmonotone_kkt_norm_ = ivalue;
    }
    else {
      nonmonotone_kkt_norm_ = 1;
    }

    if (options.GetIntegerValue("nonmonotone_kkt_centrality", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=0 && ivalue<=3, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"nonmonotone_kkt_centrality\": This value must be between 0 and 3.");
      nonmonotone_kkt_centrality_ = ivalue;
    }
    else if (options.GetIntegerValue("quality_function_centrality", ivalue, prefix)) {
      ASSERT_EXCEPTION(ivalue>=0 && ivalue<=3, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"quality_function_centrality\": This value must be between 0 and 3.");
      nonmonotone_kkt_centrality_ = ivalue;
    }
    else {
      nonmonotone_kkt_centrality_ = 0;
    }

    init_dual_inf_ = -1.;
    init_primal_inf_ = -1.;

    refs_vals_.clear();
    check_if_no_bounds_ = false;
    no_bounds_ = false;
    IpData().SetFreeMuMode(true);

    accepted_x_ = NULL;
    accepted_s_ = NULL;
    accepted_y_c_ = NULL;
    accepted_y_d_ = NULL;
    accepted_z_L_ = NULL;
    accepted_z_U_ = NULL;
    accepted_v_L_ = NULL;
    accepted_v_U_ = NULL;

    // TODO do we need to initialize the linesearch object?

    return retvalue;
  }

  void NonmonotoneMuUpdate::UpdateBarrierParameter()
  {
    // of there are not bounds, we always return the minimum MU value
    // ToDo put information on whether problem has bounds into IpCq
    if (!check_if_no_bounds_) {
      Index n_bounds = IpData().curr_z_L()->Dim() + IpData().curr_z_U()->Dim()
                       + IpData().curr_v_L()->Dim() + IpData().curr_v_U()->Dim();

      if (n_bounds==0) {
        no_bounds_ = true;
        IpData().Set_mu(mu_min_);
        IpData().Set_tau(tau_min_);
      }

      check_if_no_bounds_ = true;
    }

    if (no_bounds_)
      return;

    if (!IpData().FreeMuMode()) {
      // if we are in the fixed mu mode, we need to check if the
      // current iterate is good enough to continue with the free mode
      bool sufficient_progress = CheckSufficientProgress();
      if (sufficient_progress) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Switching back to free mu mode.\n");
        IpData().SetFreeMuMode(true);
        RememberCurrentPointAsAccepted();
      }
      else {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Remaining in fixed mu mode.\n");

        // ToDo decide whether we want this for all options
        Number sub_problem_error = IpCq().curr_barrier_error();
        Number mu = IpData().curr_mu();
        if (sub_problem_error <= kappa_epsilon_ * mu) {
          //	DBG_ASSERT(adaptive_globalization_==2);
          // If the current barrier problem has been solved sufficiently
          // well, decrease mu
          // ToDo combine this code with MonotoneMuUpdate
          Number eps_tol = IpData().epsilon_tol();

          Number new_mu = Min( kappa_mu_*mu, pow(mu, theta_mu_) );
          new_mu = Max(new_mu, eps_tol/10);
          Number new_tau = Compute_tau(mu);
          IpData().Set_mu(new_mu);
          IpData().Set_tau(new_tau);
          Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                         "Reducing mu to %e in fixed mu mode. Tau becomes %e\n", new_mu, new_tau);
          linesearch_->Reset();
        }
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
        IpData().SetFreeMuMode(false);

	if (restore_accepted_iterate_) {
	  // Restore most recent accepted iterate to start fixed mode from
	  Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
			 "Restoring most recent accepted point.\n");
	  IpData().SetTrialPrimalVariablesFromPtr(accepted_x_, accepted_s_);
	  IpData().SetTrialConstraintMultipliersFromPtr(accepted_y_c_, accepted_y_d_);
	  IpData().SetTrialBoundMultipliersFromPtr(accepted_z_L_, accepted_z_U_,
						   accepted_v_L_, accepted_v_U_);
	  IpData().AcceptTrialPoint();
	}

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

    if (IpData().FreeMuMode()) {
      // Compute the new barrier parameter via the oracle
      Number mu = free_mu_oracle_->CalculateMu();

      mu = Max(mu, mu_min_);
      Number mu_lower_safe = lower_mu_safeguard();
      if (mu < mu_lower_safe) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "mu = %e smaller than safeguard = %e. Increasing mu.\n",
                       mu, mu_lower_safe);
        mu = mu_lower_safe;
        IpData().Append_info_string("m");
      }

      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Barrier parameter mu computed by oracle is %e\n",
                     mu);

      // Apply safeguards if appropriate
      mu = Min(mu, mu_max_);
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
    if (mu_never_fix_)
      return true;

    bool retval = true;

    switch (adaptive_globalization_) {
      case 1 : {
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
      }
      break;
      case 2 : {
        retval = filter_.Acceptable(IpCq().curr_f(),
                                    IpCq().curr_constraint_violation());
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }

    return retval;
  }

  void
  NonmonotoneMuUpdate::RememberCurrentPointAsAccepted()
  {
    switch (adaptive_globalization_) {
      case 1 : {
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
      break;
      case 2 : {
        Number theta = IpCq().curr_constraint_violation();
        Number param = 1e-5;
        // ToDo need margin
        filter_.AddEntry(IpCq().curr_f() - param*theta,
                         IpCq().curr_constraint_violation() - param*theta,
                         IpData().iter_count());
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }

    if (restore_accepted_iterate_) {
      // Keep pointers to this iterate so that it could be restored
      accepted_x_ = IpData().curr_x();
      accepted_s_ = IpData().curr_s();
      accepted_y_c_ = IpData().curr_y_c();
      accepted_y_d_ = IpData().curr_y_d();
      accepted_z_L_ = IpData().curr_z_L();
      accepted_z_U_ = IpData().curr_z_U();
      accepted_v_L_ = IpData().curr_v_L();
      accepted_v_U_ = IpData().curr_v_U();
    }
  }

  Number
  NonmonotoneMuUpdate::Compute_tau(Number mu)
  {
    return Max(tau_min_, Min(1.-mu, tau_max_));
    //return tau_min_;
  }

  Number
  NonmonotoneMuUpdate::min_ref_val()
  {
    DBG_ASSERT(adaptive_globalization_==1);
    Number min_ref;
    DBG_ASSERT(refs_vals_.size()>0);
    std::list<Number>::iterator iter = refs_vals_.begin();
    min_ref = *iter;
    iter++;
    while (iter != refs_vals_.end()) {
      min_ref = Min(min_ref, *iter);
      iter++;
    }
    return min_ref;
  }

  Number
  NonmonotoneMuUpdate::max_ref_val()
  {
    DBG_ASSERT(adaptive_globalization_==1);
    Number max_ref;
    DBG_ASSERT(refs_vals_.size()>0);
    std::list<Number>::iterator iter = refs_vals_.begin();
    max_ref = *iter;
    iter++;
    while (iter != refs_vals_.end()) {
      max_ref = Max(max_ref, *iter);
      iter++;
    }
    return max_ref;
  }

  Number
  NonmonotoneMuUpdate::NewFixedMu()
  {
    Number max_ref;
    // ToDo: Decide whether we should impose an upper bound on 
    // mu based on the smallest reference value.  For now, don't
    // impose one.
    max_ref = 1e20;
    /*
    switch (adaptive_globalization_) {
      case 1 :
      max_ref = max_ref_val();
      break;
      case 2 : {
        max_ref = 1e20;
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }
    */

    Number new_mu;

    if (IsValid(fix_mu_oracle_)) {
      new_mu = fix_mu_oracle_->CalculateMu();
    }
    else {
      new_mu = fixed_mu_avrg_factor_*IpCq().curr_avrg_compl();
    }
    new_mu = Max(new_mu, lower_mu_safeguard());
    new_mu = Min(new_mu, 0.1 * max_ref);

    new_mu = Max(new_mu, mu_min_);
    new_mu = Min(new_mu, mu_max_);

    return new_mu;
  }

  //ToDo put the following into CalculatedQuantities?
  Number
  NonmonotoneMuUpdate::curr_norm_pd_system()
  {
    Index n_dual = IpData().curr_x()->Dim() + IpData().curr_s()->Dim();
    Index n_pri = IpData().curr_y_c()->Dim() + IpData().curr_y_d()->Dim();
    Index n_comp = IpData().curr_z_L()->Dim() + IpData().curr_z_U()->Dim() +
                   IpData().curr_v_L()->Dim() + IpData().curr_v_U()->Dim();

    Number dual_inf;
    Number primal_inf;
    Number complty;
    switch (nonmonotone_kkt_norm_) {
    case 1:
      dual_inf =
	IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_1);
      primal_inf =
	IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_1);
      complty =
	IpCq().curr_complementarity(0., IpoptCalculatedQuantities::NORM_1);
      dual_inf /= (Number)n_dual;
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
	primal_inf /= (Number)n_pri;
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
	complty /= (Number)n_comp;
      }
      break;
    case 2:
      dual_inf =
	IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_2);
      dual_inf *= dual_inf;
      primal_inf =
	IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_2);
      primal_inf *= primal_inf;
      complty =
	IpCq().curr_complementarity(0., IpoptCalculatedQuantities::NORM_2);
      complty *= complty;
      dual_inf /= (Number)n_dual;
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
	primal_inf /= (Number)n_pri;
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
	complty /= (Number)n_comp;
      }
      break;
    case 3:
      dual_inf =
	IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_MAX);
      primal_inf =
	IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_MAX);
      complty =
	IpCq().curr_complementarity(0., IpoptCalculatedQuantities::NORM_MAX);
      break;
    case 4:
      dual_inf =
	IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_2);
      primal_inf =
	IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_2);
      complty =
	IpCq().curr_complementarity(0., IpoptCalculatedQuantities::NORM_2);
      dual_inf /= sqrt((Number)n_dual);
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
	primal_inf /= sqrt((Number)n_pri);
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
	complty /= sqrt((Number)n_comp);
      }
      break;
    }

    Number centrality = 0.;
    if (nonmonotone_kkt_centrality_!=0) {
      Number xi = IpCq().curr_centrality_measure();
      switch (nonmonotone_kkt_centrality_) {
      case 1:
	centrality = -complty*log(xi);
	break;
      case 2:
	centrality = complty/xi;
      case 3:
	centrality = complty/pow(xi,3);
	break;
      default:
	DBG_ASSERT("Unknown value for nonmonotone_kkt_centrality_");
      }
    }
    
    Number kkt_error = primal_inf + dual_inf + complty + centrality;

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "KKT error in barrier update check:\n"
                   "  primal infeasibility: %15.6e\n"
                   "    dual infeasibility: %15.6e\n"
                   "       complementarity: %15.6e\n"
                   "            centrality: %15.6e\n"
                   "             kkt error: %15.6e\n",
                   primal_inf, dual_inf, complty, centrality, kkt_error);

    return kkt_error;
  }

  Number
  NonmonotoneMuUpdate::lower_mu_safeguard()
  {
    Number dual_inf =
      IpCq().curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_1);
    Number primal_inf =
      IpCq().curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_1);
    Index n_dual = IpData().curr_x()->Dim() + IpData().curr_s()->Dim();
    dual_inf /= (Number)n_dual;
    Index n_pri = IpData().curr_y_c()->Dim() + IpData().curr_y_d()->Dim();
    DBG_ASSERT(n_pri>0 || primal_inf==0.);
    if (n_pri>0) {
      primal_inf /= (Number)n_pri;
    }

    if (init_dual_inf_ < 0.) {
      init_dual_inf_ = Max(1., dual_inf);
    }
    if (init_primal_inf_ < 0.) {
      init_primal_inf_ = Max(1., primal_inf);
    }

    Number lower_mu_safeguard =
      Max(mu_safeguard_factor_ * (dual_inf/init_dual_inf_),
          mu_safeguard_factor_ * (primal_inf/init_primal_inf_));

    if (adaptive_globalization_==1) {
      lower_mu_safeguard = Min(lower_mu_safeguard, min_ref_val());
    }

    return lower_mu_safeguard;
  }

} // namespace Ipopt
