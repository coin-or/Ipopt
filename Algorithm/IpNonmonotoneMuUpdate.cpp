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
    Number value= 0.0;

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

    bool retvalue = mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                           options, prefix);

    // TODO do we need to initialize the linesearch object?

    return retvalue;
  }

  void NonmonotoneMuUpdate::UpdateBarrierParameter()
  {
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
    Number tau = tau_min_;
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Fraction-to-the-boundary parameter tau is %e\n",
                   tau);

    // Set the new values
    IpData().Set_mu(mu);
    IpData().Set_tau(tau);

    linesearch_->Reset();
  }

} // namespace Ipopt
