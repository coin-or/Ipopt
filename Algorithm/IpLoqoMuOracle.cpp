// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpLoqoMuOracle.hpp"

#include <limits>

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  static const Index dbg_verbosity = 0;

  LoqoMuOracle::LoqoMuOracle()
      :
      MuOracle()
  {}

  LoqoMuOracle::~LoqoMuOracle()
  {}

  bool LoqoMuOracle::InitializeImpl(const OptionsList& options,
                                    const std::string& prefix)
  {
    // The following line is only here so that
    // IpoptCalculatedQuantities::CalculateSafeSlack and the first
    // output line have something to work with
    IpData().Set_mu(1.);

    return true;
  }

  Number LoqoMuOracle::CalculateMu()
  {
    DBG_START_METH("LoqoMuOracle::CalculateMu",
                   dbg_verbosity);

    SmartPtr<const Vector> compl_x_L =  IpCq().curr_compl_x_L();
    SmartPtr<const Vector> compl_x_U =  IpCq().curr_compl_x_U();
    SmartPtr<const Vector> compl_s_L =  IpCq().curr_compl_s_L();
    SmartPtr<const Vector> compl_s_U =  IpCq().curr_compl_s_U();

    DBG_PRINT_VECTOR(2, "curr_compl_x_L", *compl_x_L);
    DBG_PRINT_VECTOR(2, "curr_compl_x_U", *compl_x_U);
    DBG_PRINT_VECTOR(2, "curr_compl_s_L", *compl_s_L);
    DBG_PRINT_VECTOR(2, "curr_compl_s_U", *compl_s_U);

    Number MinCompl = std::numeric_limits<Number>::max();
    bool have_bounds = false;

    // Compute the Minimum of all complementarities
    if( compl_x_L->Dim()>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_L->Min());
      }
      else {
        MinCompl = compl_x_L->Min();
      }
      have_bounds = true;
    }
    if( compl_x_U->Dim()>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_U->Min());
      }
      else {
        MinCompl = compl_x_U->Min();
      }
      have_bounds = true;
    }
    if( compl_s_L->Dim()>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_L->Min());
      }
      else {
        MinCompl = compl_s_L->Min();
      }
      have_bounds = true;
    }
    if( compl_s_U->Dim()>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_U->Min());
      }
      else {
        MinCompl = compl_s_U->Min();
      }
      have_bounds = true;
    }

    // If there are now bounds, just return 0.;
    if (!have_bounds) {
      return 0.;
    }

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Minimal complemantarity is %lf\n", MinCompl);

    DBG_ASSERT(MinCompl>0. && "There is a zero complementarity entry");

    Number avrg_compl = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Average complemantarity is %lf\n", avrg_compl);

    Number xi = MinCompl/avrg_compl;
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Xi (distance from uniformity) is %lf\n", xi);

    //Number factor = 1.-tau_min_;   //This is the original values
    Number factor = 0.05;   //This is the value I used otherwise
    Number mu = 0.1*pow(Min(factor*(1.-xi)/xi,2),3.)*avrg_compl;
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Barrier parameter proposed by LOQO rule is %lf\n", mu);

    // DELETEME
    char ssigma[40];
    sprintf(ssigma, "sigma=%e", 0.1*pow(Min(factor*(1.-xi)/xi,2),3.));
    IpData().Append_info_string(ssigma);

    return mu;
  }

} // namespace Ipopt
